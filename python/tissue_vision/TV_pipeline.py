#!/usr/bin/env python3

import os
import sys

from configargparse import Namespace, ArgParser
from typing import Dict, List, Union

import pandas as pd

from pydpiper.core.stages import Stages, Result
from pydpiper.core.arguments import CompoundParser, AnnotatedParser, BaseParser, _mk_lsq6_parser, to_lsq6_conf
from pydpiper.core.util import maybe_deref_path
from pydpiper.execution.application import execute, mk_application
from pydpiper.core.files import FileAtom
from pydpiper.core.arguments import application_parser, registration_parser, execution_parser, parse
from pydpiper.minc.files import MincAtom
from pydpiper.minc.registration import autocrop, create_quality_control_images, check_MINC_input_files, lsq12_nlin, \
    get_linear_configuration_from_options, LinearTransType, get_nonlinear_component, \
    get_registration_targets_from_init_model, xfmconcat, mincresample_new
from pydpiper.pipelines.MBM import mbm, mk_mbm_parser
from pydpiper.pipelines.MAGeT import Interpolation

from tissue_vision.arguments import consensus_to_atlas_parser
from tissue_vision.reconstruction import TV_stitch_wrap, cellprofiler_wrap, stacks_to_volume, \
    antsRegistration, get_like, get_through_plane_xfm, concat_xfm, mincmath

def get_imgs(options):
    if options.files:
        raise ValueError("you used --files; please use --csv-file")

    csv = pd.read_csv(options.csv_file, dtype='str', index_col=False)
    csv_base = os.path.dirname(options.csv_file)


    csv = csv.assign(anatomy_MincAtom = lambda df: df["anatomy"].apply(lambda file:
                        MincAtom(os.path.join(csv_base, file.strip()),
                                 pipeline_sub_dir=os.path.join(options.output_directory,
                                                               options.pipeline_name + "_processed"))),
                     count_MincAtom = lambda df: df["count"].apply(lambda file:
                        MincAtom(os.path.join(csv_base, file.strip()),
                                 pipeline_sub_dir=os.path.join(options.output_directory,
                                                               options.pipeline_name + "_processed"))))
    return csv

def tissue_vision_pipeline(options):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name

    csv = original_csv = get_imgs(options.application)
    # check_MINC_input_files([img.path for img in imgs])

    s = Stages()

    s.defer(create_quality_control_images(imgs=csv['anatomy_MincAtom'].tolist(), montage_dir=output_dir,
                                          montage_output=os.path.join(output_dir, pipeline_name + "_resampled",
                                                                      "input_montage"),
                                          auto_range=True,
                                          message="input_mincs"))
#############################
# Step 1: Run MBM.py to create a consensus average
#############################
    mbm_result = s.defer(mbm(imgs=csv['anatomy_MincAtom'].tolist(), options=options,
                             prefix=options.application.pipeline_name,
                             output_dir=output_dir,
                             with_maget=False))

    #TODO remove
    transforms = mbm_result.xfms.assign(native_file=lambda df: df.rigid_xfm.apply(lambda x: x.source),
                            lsq6_file=lambda df: df.lsq12_nlin_xfm.apply(lambda x: x.source),
                            lsq6_mask_file=lambda df:
                              df.lsq12_nlin_xfm.apply(lambda x: x.source.mask if x.source.mask else ""),
                            nlin_file=lambda df: df.lsq12_nlin_xfm.apply(lambda x: x.resampled),
                            nlin_mask_file=lambda df:
                              df.lsq12_nlin_xfm.apply(lambda x: x.resampled.mask if x.resampled.mask else ""))\
        .applymap(maybe_deref_path)
    determinants = mbm_result.determinants.drop(["full_det", "nlin_det"], axis=1)\
        .applymap(maybe_deref_path)

    csv = csv.assign(anatomy_lsq6_MincAtom=mbm_result.xfms.lsq12_nlin_xfm.apply(lambda xfm: xfm.source),
                     mbm_lsq6_XfmAtom=mbm_result.xfms.rigid_xfm.apply(lambda x: x.xfm),
                     mbm_lsq12_nlin_XfmAtom=mbm_result.xfms.lsq12_nlin_xfm.apply(lambda x: x.xfm),
                     mbm_full_XfmAtom=mbm_result.xfms.overall_xfm.apply(lambda x: x.xfm))

    # x.assign(count_lsq6_MincAtom=lambda df: [x + y for x, y in zip(df["x"], df["y"])])
    csv = csv.assign(count_lsq6_MincAtom = lambda df:
    [s.defer(mincresample_new(img = img,
                              xfm = xfm,
                              like = like))
     for img, xfm, like in zip(df["count_MincAtom"],
                               df["mbm_lsq6_XfmAtom"],
                               df["anatomy_lsq6_MincAtom"])])


#############################
# Step 2: Register consensus average to ABI tissuevision Atlas
#############################
    lsq12_conf = get_linear_configuration_from_options(conf=options.mbm.lsq12,
                                                       transform_type=LinearTransType.lsq12,
                                                       file_resolution=options.registration.resolution)
    nlin_component = get_nonlinear_component(reg_method=options.mbm.nlin.reg_method)

    atlas_target = MincAtom(name=options.consensus_to_atlas.atlas_target,
                            orig_name=options.consensus_to_atlas.atlas_target,
                            mask=MincAtom(name=options.consensus_to_atlas.atlas_target_mask,
                                          orig_name=options.consensus_to_atlas.atlas_target_mask))
    atlas_target_label = MincAtom(name=options.consensus_to_atlas.atlas_target_label,
                                  orig_name=options.consensus_to_atlas.atlas_target_label,
                                  mask=MincAtom(name=options.consensus_to_atlas.atlas_target_mask,
                                                orig_name=options.consensus_to_atlas.atlas_target_mask))

    lsq12_nlin_result = s.defer(lsq12_nlin(source=mbm_result.avg_img,
                                           target=atlas_target,
                                           lsq12_conf=lsq12_conf,
                                           nlin_module=nlin_component,
                                           nlin_options=options.mbm.nlin.nlin_protocol,
                                           resolution=options.registration.resolution,
                                           resample_source=False
                                           ))

#############################
# Step 3: Resample count volumes to ABI tissuevision Atlas space and vice versa
#############################

    csv = csv.assign(lsq6_to_atlas_XfmAtom = lambda df: df['mbm_lsq12_nlin_XfmAtom'].apply(lambda xfm:
                            s.defer(xfmconcat([xfm, lsq12_nlin_result.xfm]))))

    csv = csv.assign(
        anatomy_targetspace_MincAtom=lambda df:
        [s.defer(mincresample_new(img=img, xfm=xfm, like=atlas_target))
         for img, xfm in zip(df["anatomy_lsq6_MincAtom"], df["lsq6_to_atlas_XfmAtom"])],
        count_targetspace_MincAtom=lambda df:
        [s.defer(mincresample_new(img=img, xfm=xfm, like=atlas_target))
         for img, xfm in zip(df["count_lsq6_MincAtom"], df["lsq6_to_atlas_XfmAtom"])],
        atlas_lsq6space_MincAtom=lambda df:
        [s.defer(mincresample_new(img=atlas_target_label, xfm=xfm, like=like, invert=True,
                                  interpolation=Interpolation.nearest_neighbour,
                                  extra_flags=('-keep_real_range',)))
         for xfm, like in zip( df["lsq6_to_atlas_XfmAtom"], df["count_lsq6_MincAtom"])]
    )

    csv.applymap(maybe_deref_path).to_csv("analysis.csv",index=False)

    s.defer(create_quality_control_images(imgs=csv.count_targetspace_MincAtom.tolist(), montage_dir=output_dir,
                                          montage_output=os.path.join(output_dir, pipeline_name + "_resampled",
                                                                      "count_montage"),
                                          auto_range=True,
                                          message="count_mincs"))
    return Result(stages=s, output=())


lsq6_parser = AnnotatedParser(parser=BaseParser(_mk_lsq6_parser(with_nuc=False, with_inormalize=False), "LSQ6"),
                              namespace="lsq6", cast=to_lsq6_conf)
tv_pipeline_application = mk_application(parsers=[AnnotatedParser(parser=mk_mbm_parser(with_common_space=False,
                                                                                       lsq6_parser=lsq6_parser),
                                                                  namespace="mbm"),
                                                  consensus_to_atlas_parser],
                                         pipeline=tissue_vision_pipeline)

if __name__ == "__main__":
    tv_pipeline_application()
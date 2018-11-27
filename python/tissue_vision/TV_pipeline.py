#!/usr/bin/env python3

import os
import sys

from configargparse import Namespace, ArgParser
from typing import Dict, List, Union

import pandas as pd

from pydpiper.core.stages import Stages, Result
from pydpiper.core.arguments import CompoundParser, AnnotatedParser, BaseParser
from pydpiper.core.util import maybe_deref_path
from pydpiper.execution.application import execute, mk_application
from pydpiper.core.files import FileAtom
from pydpiper.core.arguments import application_parser, registration_parser, execution_parser, parse
from pydpiper.minc.files import MincAtom
from pydpiper.minc.registration import autocrop, create_quality_control_images, check_MINC_input_files, lsq12_nlin, \
    get_linear_configuration_from_options, LinearTransType, get_nonlinear_component, \
    get_registration_targets_from_init_model, xfmconcat
from pydpiper.pipelines.MBM import mbm, mk_mbm_parser

from tissue_vision.arguments import TV_stitch_parser, cellprofiler_parser, stacks_to_volume_parser, autocrop_parser

from tissue_vision.reconstruction import TV_stitch_wrap, cellprofiler_wrap, stacks_to_volume, \
    antsRegistration, get_like, tif_to_minc, get_through_plane_xfm, concat_xfm, mincresample, mincmath
from tissue_vision.TV_stitch import get_params

def tissue_vision_pipeline(options):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name

    s = Stages()

#############################
# Step 1: Run MBM.py to create a consensus average
#############################
    #this is annoying but consistent with parser logic
    if "--run-mbm" in sys.argv[1:]:
        imgs = all_smooth_pad_results

        check_MINC_input_files([img.path for img in imgs])

        mbm_result = s.defer(mbm(imgs=imgs, options=options,
                                 prefix=options.application.pipeline_name,
                                 output_dir=output_dir,
                                 with_maget=False))

        # create useful CSVs (note the files listed therein won't yet exist ...):
        transforms = mbm_result.xfms.assign(native_file=lambda df: df.rigid_xfm.apply(lambda x: x.source),
                                lsq6_file=lambda df: df.lsq12_nlin_xfm.apply(lambda x: x.source),
                                lsq6_mask_file=lambda df:
                                  df.lsq12_nlin_xfm.apply(lambda x: x.source.mask if x.source.mask else ""),
                                nlin_file=lambda df: df.lsq12_nlin_xfm.apply(lambda x: x.resampled),
                                nlin_mask_file=lambda df:
                                  df.lsq12_nlin_xfm.apply(lambda x: x.resampled.mask if x.resampled.mask else ""))\
            .applymap(maybe_deref_path)
        transforms.to_csv("transforms.csv", index=False)

        determinants = mbm_result.determinants.drop(["full_det", "nlin_det"], axis=1)\
            .applymap(maybe_deref_path)
        determinants.to_csv("determinants.csv", index=False)



#############################
# Step 2: Register consensus average to ABI tissuevision Atlas
#############################
        lsq12_conf = get_linear_configuration_from_options(conf=options.mbm.lsq12,
                                                           transform_type=LinearTransType.lsq12,
                                                           file_resolution=options.registration.resolution)
        nlin_component = get_nonlinear_component(reg_method=options.mbm.nlin.reg_method)

        atlas_target = MincAtom(name=options.tissue_vision.final_registration.atlas_target,
                                orig_name=options.tissue_vision.final_registration.atlas_target,
                                mask=MincAtom(name=options.tissue_vision.final_registration.atlas_target_mask,
                                              orig_name=options.tissue_vision.final_registration.atlas_target_mask))

        lsq12_nlin_result = s.defer(lsq12_nlin(source=mbm_result.avg_img,
                                               target=atlas_target,
                                               lsq12_conf=lsq12_conf,
                                               nlin_module=nlin_component,
                                               nlin_options=options.mbm.nlin.nlin_protocol,
                                               resolution=options.registration.resolution,
                                               resample_source=False
                                               ))

#############################
# Step 3: Resample binary volumes to ABI tissuevision Atlas space and vice versa
#############################
        all_full_xfms = []
        init_model = get_registration_targets_from_init_model(init_model_standard_file=options.mbm.lsq6.target_file,
                                                              output_dir=output_dir,
                                                              pipeline_name=pipeline_name)
        all_binary_pad_lsq6 = [xfm_handler._resampled for xfm_handler in mbm_result.xfms.rigid_xfm]
        for mbm_lsq12_nlin_xfm, binary_pad_lsq6, binary_resampled, atlas_resampled in \
                zip(mbm_result.xfms.lsq12_nlin_xfm, all_binary_pad_lsq6, all_binary_resampled, all_atlas_resampled):
            full_xfm = s.defer(xfmconcat([mbm_lsq12_nlin_xfm.xfm, lsq12_nlin_result.xfm]))
            all_full_xfms.append(full_xfm)

            s.defer(mincresample(img=binary_pad_lsq6,
                                 xfm=full_xfm,
                                 like=atlas_target,
                                 resampled=binary_resampled,
                                 output_dir=output_dir))

            s.defer(mincresample(img=atlas_target,
                                 xfm=full_xfm,
                                 like=binary_pad_lsq6,
                                 resampled=atlas_resampled,
                                 output_dir=output_dir,
                                 invert_xfm=True))

        transforms=transforms.assign(binary_resampled=[minc_atom.path for minc_atom in all_binary_resampled],
                                     atlas_resampled = [minc_atom.path for minc_atom in all_atlas_resampled],
                                     full_xfm=[xfm.path for xfm in all_full_xfms])

        analysis = transforms.merge(determinants, left_on="lsq12_nlin_xfm", right_on="inv_xfm", how='inner') \
            .drop(["xfm", "inv_xfm"], axis=1)
        reconstructed.merge(analysis, left_on="smooth_padded", right_on="native_file").drop(["native_file"], axis=1)\
            .to_csv("analysis.csv",index=False)

        s.defer(create_quality_control_images(imgs=all_binary_resampled, montage_dir=output_dir,
                                              montage_output=os.path.join(output_dir, pipeline_name + "_resampled",
                                                                          "%s_montage" % binary),
                                              auto_range=True,
                                              message="%s_mincs" % binary))
    return Result(stages=s, output=())

#############################
# Combine Parser & Make Application
#############################
def mk_tissue_vision_parser():
    p = ArgParser(add_help=False)
    p.add_argument("--atlas-target", dest="atlas_target",
                   type=str,
                   default=None,
                   help="Register the consensus average to the ABI Atlas")
    p.add_argument("--atlas-target-mask", dest="atlas_target_mask",
                   type=str,
                   default=None,
                   help="Register the consensus average to the ABI Atlas")
    registration_parser = AnnotatedParser(parser=BaseParser(p, 'tissue_vision'), namespace='final_registration')
    p = ArgParser(add_help=False)
    p.add_argument("--run-mbm", dest="run_mbm",
                   action="store_true", default=False,
                   help="Run MBM after reconstructing the brains. Use this flag with --help to get MBM help.")
    mbm_parser = AnnotatedParser(parser=BaseParser(p, 'mbm'), namespace='mbm')
    return CompoundParser([TV_stitch_parser,
              cellprofiler_parser,
              stacks_to_volume_parser,
              autocrop_parser, registration_parser, mbm_parser])

#############################
# Run
#############################

if __name__ == "__main__":
    #this is needed as adding the parser without passing the flags breaks
    if "--run-mbm" in sys.argv[1:]:
        p = CompoundParser([application_parser, registration_parser, execution_parser,
                            AnnotatedParser(parser=mk_tissue_vision_parser(), namespace='tissue_vision'),
                            AnnotatedParser(parser=mk_mbm_parser(with_common_space=False, with_maget=False),
                                            namespace="mbm")])

        # hacky index-based reaching into the mbm-lsq6 parser to turn off
        p.parsers[-1].parser.parsers[0].parser.argparser.set_defaults(inormalize=False)
        p.parsers[-1].parser.parsers[0].parser.argparser.set_defaults(nuc=False)

    else:
        p = CompoundParser([application_parser, registration_parser, execution_parser,
                           AnnotatedParser(parser=mk_tissue_vision_parser(), namespace='tissue_vision')])

    parsed_options = parse(p, sys.argv[1:])
    execute(tissue_vision_pipeline(parsed_options).stages, parsed_options)

lsq6_application = mk_application(parsers=[lsq6_parser], pipeline=lsq6_pipeline)

if __name__ == "__main__":
    lsq6_application()
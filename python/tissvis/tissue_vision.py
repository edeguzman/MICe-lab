#!/usr/bin/env python3

import os

from configargparse import Namespace
from typing import Dict, List, Union

import pandas as pd

from pydpiper.core.stages import Stages, Result
from pydpiper.core.arguments import CompoundParser, AnnotatedParser
from pydpiper.core.util import maybe_deref_path
from pydpiper.execution.application import mk_application
from pydpiper.core.files import FileAtom
from pydpiper.minc.files import MincAtom
from pydpiper.minc.registration import autocrop, check_MINC_input_files
from pydpiper.pipelines.MBM import mbm, MBMConf, common_space, mk_mbm_parser

from tissvis.arguments import TV_stitch_parser, cellprofiler_parser, stacks_to_volume_parser, autocrop_parser
from tissvis.reconstruction import TV_stitch_wrap, cellprofiler_wrap, stacks_to_volume, \
    antsRegistration, get_like, tif_to_minc, get_through_plane_xfm, concat_xfm, mincresample, mincmath
from tissvis.TV_stitch import get_params

class Brain(object):
    def __init__(self,
                 brain_directory: FileAtom,
                 name: str,
                 z_start: Union[int,None],
                 z_end: Union[int,None],
                 z_section: Union[int,None],
                 ) -> None:
        self.brain_directory = brain_directory
        self.name = name
        self.z_start = z_start
        self.z_end = z_end
        self.z_section = z_section

def get_brains(options):
    if options.files:
        raise ValueError("you used --files; please use --csv-file")

    csv = pd.read_csv(options.csv_file)

    brains = [Brain(FileAtom(brain_directory), brain_name, z_start, z_end, z_section)
              for brain_directory, brain_name, z_start, z_end,z_section
              in zip(csv.brain_directory, csv.brain_name, csv.Zstart, csv.Zend, csv.Zsection)]
    return brains


def tissue_vision_pipeline(options):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name

    s = Stages()

    brains = get_brains(options.application) # List(Brain,...)

    env_vars = {}
    env_vars['PYTHONPATH'] = options.tissue_vision.cellprofiler.python2_path
    cppline = FileAtom(options.tissue_vision.cellprofiler.cellprofiler_pipeline)

    # Hold results obtained in the loop
    all_TV_stitch_results = []
    all_cellprofiler_results = []
    all_binary_volume_results = []
    all_smooth_volume_results = []
    all_binary_volume_isotropic_results = []
    all_smooth_volume_isotropic_results = []
    all_smooth_pad_results = []
    all_binary_pad_results = []

#############################
# Step 1: Run TV_stitch.py
#############################
    for brain in brains:
        slice_directory = os.path.join(output_dir, pipeline_name + "_stitched", brain.name)

        stitched = []
        brain.x, brain.y, brain.z, brain.z_resolution = get_params(os.path.join(brain.brain_directory.path, brain.name))

        #accounting for errors in tile acquisition
        brain.z_start = 1 if pd.isna(brain.z_start) else int(brain.z_start)
        brain.z_end = brain.z if pd.isna(brain.z_end) else int(brain.z_end)
        brain.z_section = None if pd.isna(brain.z_section) else int(brain.z_section)

        for z in range (brain.z_start, brain.z_end + 1):
            brain.slice_stitched = FileAtom(os.path.join(slice_directory, brain.name + "_Z%04d.tif" % z))
            stitched.append(brain.slice_stitched)

        if not brain.z_section:
            TV_stitch_result = s.defer(TV_stitch_wrap(brain_directory = brain.brain_directory,
                                                    brain_name = brain.name,
                                                    stitched = stitched,
                                                    TV_stitch_options = options.tissue_vision.TV_stitch,
                                                    Zstart=brain.z_start,
                                                    Zend=brain.z_end,
                                                    output_dir = output_dir
                                                    ))
            all_TV_stitch_results.append(TV_stitch_result)

        if brain.z_section:
            TV_stitch_result = s.defer(TV_stitch_wrap(brain_directory=brain.brain_directory,
                                                      brain_name=brain.name,
                                                      stitched=stitched[brain.z_start:brain.z_section],
                                                      TV_stitch_options=options.tissue_vision.TV_stitch,
                                                      Zstart=brain.z_start,
                                                      Zend=brain.z_section - 1,
                                                      output_dir=output_dir
                                                      ))
            all_TV_stitch_results.append(TV_stitch_result)

            TV_stitch_result = s.defer(TV_stitch_wrap(brain_directory=brain.brain_directory,
                                                      brain_name=brain.name,
                                                      stitched=stitched[brain.z_section:brain.z_end-1],
                                                      TV_stitch_options=options.tissue_vision.TV_stitch,
                                                      Zstart=brain.z_section,
                                                      Zend=brain.z_end,
                                                      output_dir=output_dir
                                                      ))
            all_TV_stitch_results.append(TV_stitch_result)

#TODO write a when_finished_hook to tell the user that this finished.
#############################
# Step 2: Run cellprofiler
#############################
        anatomical = options.tissue_vision.cellprofiler.anatomical_name
        binary = options.tissue_vision.cellprofiler.binary_name

        brain.cp_directory = os.path.join(output_dir, pipeline_name + "_cellprofiler", brain.name)
        brain.batch_data = FileAtom(os.path.join(brain.cp_directory,"Batch_data.h5"))

        overLays = []
        smooths = []
        binaries = []

        for z in range(brain.z_start, brain.z_end+1):
            brain.slice_overLay = FileAtom(
                os.path.join(brain.cp_directory, brain.name + "_Z%04d_overLay.tiff" % z))
            brain.slice_smooth = FileAtom(
                os.path.join(brain.cp_directory, brain.name + "_Z%04d_" % z + anatomical + ".tiff"))
            brain.slice_binary = FileAtom(
                os.path.join(brain.cp_directory, brain.name + "_Z%04d_" % z + binary + ".tiff"))
            overLays.append(brain.slice_overLay)
            smooths.append(brain.slice_smooth)
            binaries.append(brain.slice_binary)

        cellprofiler_result = s.defer(cellprofiler_wrap(stitched = stitched,
                                                              cellprofiler_pipeline = cppline,
                                                              batch_data = brain.batch_data,
                                                              overLays = overLays,
                                                              smooths = smooths,
                                                              binaries = binaries,
                                                              Zstart = brain.z_start,
                                                              Zend = brain.z_end,
                                                              output_dir = output_dir,
                                                              env_vars = env_vars
                                                             ))
        all_cellprofiler_results.append(cellprofiler_result)

#############################
# Step 3: Run stacks_to_volume.py
#############################
        smooth_volume = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + anatomical + "_stacked.mnc"))
        binary_volume = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                      brain.name + "_" + binary + "_stacked.mnc"))

        if not brain.z_section:
            smooth_slices_to_volume_results = s.defer(stacks_to_volume(
                slices = smooths,
                volume = smooth_volume,
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                uniform_sum=False,
                z_resolution=brain.z_resolution,
                output_dir=output_dir
                ))
            all_smooth_volume_results.append(smooth_slices_to_volume_results)

            binary_slices_to_volume_results = s.defer(stacks_to_volume(
                slices = binaries,
                volume = binary_volume,
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                z_resolution=brain.z_resolution,
                uniform_sum = True,
                output_dir=output_dir
                ))
            all_binary_volume_results.append(binary_slices_to_volume_results)

        if brain.z_section:

            smooth_volume_1 = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + anatomical + "_stacked_1.mnc"))
            smooth_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=smooths[brain.z_start-1:brain.z_section-1],
                volume=smooth_volume_1,
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                uniform_sum=False,
                z_resolution=brain.z_resolution,
                output_dir=output_dir
            ))
            all_smooth_volume_results.append(smooth_slices_to_volume_results)

            smooth_volume_2 = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + anatomical + "_stacked_2.mnc"))
            smooth_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=smooths[brain.z_section-1:brain.z_end],
                volume=smooth_volume_2,
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                uniform_sum=False,
                z_resolution=brain.z_resolution,
                output_dir=output_dir
            ))
            all_smooth_volume_results.append(smooth_slices_to_volume_results)


            binary_volume_1 = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + binary + "_stacked_1.mnc"))
            binary_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=binaries[brain.z_start-1:brain.z_section-1],
                volume=binary_volume_1,
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                z_resolution=brain.z_resolution,
                uniform_sum=True,
                output_dir=output_dir
            ))
            all_binary_volume_results.append(binary_slices_to_volume_results)

            binary_volume_2 = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + binary + "_stacked_2.mnc"))
            binary_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=binaries[brain.z_section-1:brain.z_end],
                volume=binary_volume_2,
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                z_resolution=brain.z_resolution,
                uniform_sum=True,
                output_dir=output_dir
            ))
            all_binary_volume_results.append(binary_slices_to_volume_results)

            #create minc slices for registration
            target_minc = s.defer(tif_to_minc(
                tif=smooths[brain.z_section - 2],
                volume=MincAtom(smooths[brain.z_section - 2].path.replace("tiff","mnc")),
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                z_resolution=brain.z_resolution,
                output_dir=output_dir))
            source_minc = s.defer(tif_to_minc(
                tif=smooths[brain.z_section - 1],
                volume=MincAtom(smooths[brain.z_section - 1].path.replace("tiff","mnc")),
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                z_resolution=brain.z_resolution,
                output_dir=output_dir))

            #create in-plane transform
            in_plane_transform = FileAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "0_GenericAffine.xfm"))
            s.defer(antsRegistration(img = source_minc,
                                     target = target_minc,
                                     transform = in_plane_transform,
                                     output_dir=output_dir))

            #create through-plane transform
            through_plane_xfm = FileAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_through_plane.xfm"))
            s.defer(get_through_plane_xfm(img = smooth_volume_1,
                                          xfm = through_plane_xfm,
                                          output_dir = output_dir))

            #concatenate the transforms
            xfm_concat = FileAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_concat.xfm"))
            s.defer(concat_xfm(xfms=[in_plane_transform, through_plane_xfm],
                               outxfm = xfm_concat,
                               output_dir = output_dir))

            # create like file from first section
            smooth_volume_1_like = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                         brain.name + "_" + anatomical + "_stacked_like.mnc"))
            s.defer(get_like(img=smooth_volume_1, ref=smooth_volume_2,
                             like=smooth_volume_1_like, output_dir=output_dir))
            binary_volume_1_like = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                         brain.name + "_" + binary + "_stacked_like.mnc"))
            s.defer(get_like(img=binary_volume_1, ref=binary_volume_2,
                             like=binary_volume_1_like, output_dir=output_dir))

            #transform the second section
            smooth_volume_2_transformed = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                brain.name + "_" + anatomical + "_stacked_2_transformed.mnc"))
            s.defer(mincresample(img = smooth_volume_2,
                                 xfm = xfm_concat,
                                 like = smooth_volume_1_like,
                                 resampled = smooth_volume_2_transformed,
                                 output_dir = output_dir))
            binary_volume_2_transformed = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                    brain.name + "_" + binary + "_stacked_2_transformed.mnc"))
            s.defer(mincresample(img=binary_volume_2,
                                 xfm=xfm_concat,
                                 like=binary_volume_1_like,
                                 resampled=binary_volume_2_transformed,
                                 output_dir=output_dir))

            #add the first section's like with the transformed second section
            s.defer(mincmath(imgs = [smooth_volume_1_like, smooth_volume_2_transformed],
                             result = smooth_volume,
                             output_dir=output_dir))
            s.defer(mincmath(imgs=[binary_volume_1_like, binary_volume_2_transformed],
                             result=binary_volume,
                             output_dir=output_dir))

#############################
# Step 4: Run autocrop to resample to isotropic
#############################
        smooth_volume_isotropic = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                        brain.name + "_" + anatomical + "_stacked_isotropic.mnc"))
        smooth_volume_isotropic_results = s.defer(autocrop(
            isostep = options.tissue_vision.stacks_to_volume.plane_resolution,
            img = smooth_volume,
            autocropped = smooth_volume_isotropic
        ))
        all_smooth_volume_isotropic_results.append(smooth_volume_isotropic_results)

        binary_volume_isotropic = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                 brain.name + "_" + binary + "_stacked_isotropic.mnc"))
        binary_volume_isotropic_results = s.defer(autocrop(
            isostep = options.tissue_vision.stacks_to_volume.plane_resolution,
            img = binary_volume,
            autocropped = binary_volume_isotropic,
            nearest_neighbour = True
        ))
        all_binary_volume_isotropic_results.append(binary_volume_isotropic_results)

#############################
# Step 5: Run autocrop to pad the isotropic images
#############################
        x_pad = options.tissue_vision.autocrop.x_pad
        y_pad = options.tissue_vision.autocrop.y_pad
        z_pad = options.tissue_vision.autocrop.z_pad

        smooth_padded = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                        brain.name + "_" + anatomical + "_padded.mnc"))
        smooth_pad_results = s.defer(autocrop(
            img = smooth_volume_isotropic,
            autocropped = smooth_padded,
            x_pad = x_pad,
            y_pad = y_pad,
            z_pad = z_pad
        ))
        all_smooth_pad_results.append(smooth_pad_results)

        binary_padded = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                        brain.name + "_" + binary + "_padded.mnc"))
        binary_pad_results = s.defer(autocrop(
            img = binary_volume_isotropic,
            autocropped = binary_padded,
            x_pad = x_pad,
            y_pad = y_pad,
            z_pad = z_pad
        ))
        all_binary_pad_results.append(binary_pad_results)

#############################
# Step 6: Run MBM.py
#############################
    # #TODO turn off inormalize and nuc from the parser...somehow
    # check_MINC_input_files([img.path for img in all_smooth_pad_results])
    #
    # mbm_result = s.defer(mbm(imgs=all_smooth_pad_results, options=options,
    #                          prefix=pipeline_name, output_dir=output_dir))
    #
    # if options.mbm.common_space.do_common_space_registration:
    #     s.defer(common_space(mbm_result, options))
    #
    # # create useful CSVs (note the files listed therein won't yet exist ...):
    # #TODO why is this so tediously the same as mbm_pipeline()
    # (mbm_result.xfms.assign(native_file=lambda df: df.rigid_xfm.apply(lambda x: x.source),
    #                         lsq6_file=lambda df: df.lsq12_nlin_xfm.apply(lambda x: x.source),
    #                         lsq6_mask_file=lambda df:
    #                           df.lsq12_nlin_xfm.apply(lambda x: x.source.mask if x.source.mask else ""),
    #                         nlin_file=lambda df: df.lsq12_nlin_xfm.apply(lambda x: x.resampled),
    #                         common_space_file=lambda df: df.xfm_to_common.apply(lambda x: x.resampled)
    #                                             if options.mbm.common_space.do_common_space_registration else None)
    #  .applymap(maybe_deref_path)
    #  .drop(["common_space_file"] if not options.mbm.common_space.do_common_space_registration else [], axis=1)
    #  .to_csv("transforms.csv", index=False))
    #
    # (mbm_result.determinants.drop(["full_det", "nlin_det"], axis=1)
    #  .applymap(maybe_deref_path).to_csv("determinants.csv", index=False))

    return Result(stages=s, output=())

#############################
# Combine Parser & Make Application
#############################
def mk_tissue_vision_parser():
    return CompoundParser([TV_stitch_parser,
              cellprofiler_parser,
              stacks_to_volume_parser,
              autocrop_parser])

#tissue_vision_parser = AnnotatedParser(CompoundParser([TV_stitch_parser, cellprofiler_parser, stacks_to_volume_parser,
                                       # autocrop_parser]), namespace='tissue_vision')

tissue_vision_application = mk_application(
    parsers=[AnnotatedParser(parser=mk_tissue_vision_parser(), namespace='tissue_vision')],
    pipeline=tissue_vision_pipeline)

#############################
# Run
#############################

if __name__ == "__main__":
    tissue_vision_application()

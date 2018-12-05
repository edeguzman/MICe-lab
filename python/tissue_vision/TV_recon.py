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

    csv = pd.read_csv(options.csv_file, dtype='str')

    brains = [Brain(FileAtom(brain_directory), brain_name, z_start, z_end, z_section)
              for brain_directory, brain_name, z_start, z_end,z_section
              in zip(csv.brain_directory, csv.brain_name, csv.Zstart, csv.Zend, csv.Zsection)]
    return brains


def tv_recon_pipeline(options):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name

    s = Stages()

    brains = get_brains(options.application) # List(Brain,...)

    cppline = FileAtom(options.cellprofiler.cellprofiler_pipeline)

    # Hold results obtained in the loop
    all_anatomical_pad_results = []
    all_count_pad_results = []
    reconstructed_mincs = []
    all_count_resampled = []
    all_atlas_resampled = []

#############################
# Step 1: Run TV_stitch.py
#############################
    for brain in brains:
        slice_directory = os.path.join(output_dir, pipeline_name + "_stitched", brain.name)

        stitched = []
        brain.x, brain.y, brain.z, brain.z_resolution = \
            get_params(os.path.join(brain.brain_directory.path, brain.name))

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
                                                    TV_stitch_options = options.TV_stitch,
                                                    Zstart=brain.z_start,
                                                    Zend=brain.z_end,
                                                    output_dir = output_dir
                                                    ))

        if brain.z_section:
            TV_stitch_result = s.defer(TV_stitch_wrap(brain_directory=brain.brain_directory,
                                                      brain_name=brain.name,
                                                      stitched=stitched[0 : brain.z_section - brain.z_start],
                                                      TV_stitch_options=options.TV_stitch,
                                                      Zstart=brain.z_start,
                                                      Zend=brain.z_section - 1,
                                                      output_dir=output_dir
                                                      ))

            TV_stitch_result = s.defer(TV_stitch_wrap(brain_directory=brain.brain_directory,
                                                      brain_name=brain.name,
                                                      stitched=stitched[brain.z_section - brain.z_start:brain.z_end-1],
                                                      TV_stitch_options=options.TV_stitch,
                                                      Zstart=brain.z_section,
                                                      Zend=brain.z_end,
                                                      output_dir=output_dir
                                                      ))

#TODO write a when_finished_hook to tell the user that this finished.
#############################
# Step 2: Run cellprofiler
#############################
        anatomical = options.cellprofiler.anatomical_name
        count = options.cellprofiler.count_name

        brain.cp_directory = os.path.join(output_dir, pipeline_name + "_cellprofiler", brain.name)
        brain.batch_data = FileAtom(os.path.join(brain.cp_directory,"Batch_data.h5"))

        overLays = []
        anatomicals = []
        counts = []

        for z in range(brain.z_start, brain.z_end+1):
            brain.slice_overLay = FileAtom(
                os.path.join(brain.cp_directory, brain.name + "_Z%04d_overLay.tiff" % z))
            brain.slice_anatomical = FileAtom(
                os.path.join(brain.cp_directory, brain.name + "_Z%04d_" % z + anatomical + ".tiff"))
            brain.slice_count = FileAtom(
                os.path.join(brain.cp_directory, brain.name + "_Z%04d_" % z + count + ".tiff"))
            overLays.append(brain.slice_overLay)
            anatomicals.append(brain.slice_anatomical)
            counts.append(brain.slice_count)

        cellprofiler_result = s.defer(cellprofiler_wrap(stitched = stitched,
                                                              cellprofiler_pipeline = cppline,
                                                              batch_data = brain.batch_data,
                                                              overLays = overLays,
                                                              anatomicals = anatomicals,
                                                              count = counts,
                                                              Zstart = brain.z_start,
                                                              Zend = brain.z_end,
                                                              output_dir = output_dir
                                                             ))

#############################
# Step 3: Run stacks_to_volume.py
#############################
        anatomical_volume = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + anatomical + "_stacked.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
        count_volume = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                      brain.name + "_" + count + "_stacked.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))

        if not brain.z_section:
            anatomical_slices_to_volume_results = s.defer(stacks_to_volume(
                slices = anatomicals,
                volume = anatomical_volume,
                stacks_to_volume_options=options.stacks_to_volume,
                uniform_sum=False,
                z_resolution=brain.z_resolution,
                output_dir=output_dir
                ))

            count_slices_to_volume_results = s.defer(stacks_to_volume(
                slices = counts,
                volume = count_volume,
                stacks_to_volume_options=options.stacks_to_volume,
                z_resolution=brain.z_resolution,
                uniform_sum = True,
                output_dir=output_dir
                ))

        if brain.z_section:

            anatomical_volume_1 = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + anatomical + "_stacked_1.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            anatomical_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=anatomicals[0 : brain.z_section - brain.z_start],
                volume=anatomical_volume_1,
                stacks_to_volume_options=options.stacks_to_volume,
                uniform_sum=False,
                z_resolution=brain.z_resolution,
                output_dir=output_dir
            ))

            anatomical_volume_2 = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + anatomical + "_stacked_2.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            anatomical_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=anatomicals[brain.z_section - brain.z_start:brain.z_end-1],
                volume=anatomical_volume_2,
                stacks_to_volume_options=options.stacks_to_volume,
                uniform_sum=False,
                z_resolution=brain.z_resolution,
                output_dir=output_dir
            ))

            count_volume_1 = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + count + "_stacked_1.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            count_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=counts[0 : brain.z_section - brain.z_start],
                volume=count_volume_1,
                stacks_to_volume_options=options.stacks_to_volume,
                z_resolution=brain.z_resolution,
                uniform_sum=True,
                output_dir=output_dir
            ))

            count_volume_2 = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + count + "_stacked_2.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            count_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=counts[brain.z_section - brain.z_start:brain.z_end-1],
                volume=count_volume_2,
                stacks_to_volume_options=options.stacks_to_volume,
                z_resolution=brain.z_resolution,
                uniform_sum=True,
                output_dir=output_dir
            ))

            #create minc slices for registration
            fixed_minc = s.defer(tif_to_minc(
                tif=anatomicals[brain.z_section - brain.z_start - 1],
                volume=MincAtom(anatomicals[brain.z_section - brain.z_start - 1].path.replace("tiff","mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked")),
                stacks_to_volume_options=options.stacks_to_volume,
                z_resolution=brain.z_resolution,
                output_dir=output_dir))
            moving_minc = s.defer(tif_to_minc(
                tif=anatomicals[brain.z_section - brain.z_start],
                volume=MincAtom(anatomicals[brain.z_section - brain.z_start].path.replace("tiff","mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked")),
                stacks_to_volume_options=options.stacks_to_volume,
                z_resolution=brain.z_resolution,
                output_dir=output_dir))

            #create in-plane transform
            in_plane_transform = FileAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "0_GenericAffine.xfm"))
            s.defer(antsRegistration(fixed = fixed_minc,
                                     moving = moving_minc,
                                     transform = in_plane_transform,
                                     output_dir=output_dir))

            #create through-plane transform
            through_plane_xfm = FileAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_through_plane.xfm"))
            s.defer(get_through_plane_xfm(img = anatomical_volume_1,
                                          xfm = through_plane_xfm,
                                          output_dir = output_dir))

            #concatenate the transforms
            xfm_concat = FileAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_concat.xfm"))
            s.defer(concat_xfm(xfms=[in_plane_transform, through_plane_xfm],
                               outxfm = xfm_concat,
                               output_dir = output_dir))

            # create like file from first section
            anatomical_volume_1_like = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                         brain.name + "_" + anatomical + "_stacked_like.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            s.defer(get_like(img=anatomical_volume_1, ref=anatomical_volume_2,
                             like=anatomical_volume_1_like, output_dir=output_dir))
            count_volume_1_like = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                         brain.name + "_" + count + "_stacked_like.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            s.defer(get_like(img=count_volume_1, ref=count_volume_2,
                             like=count_volume_1_like, output_dir=output_dir))

            #transform the second section
            anatomical_volume_2_transformed = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                brain.name + "_" + anatomical + "_stacked_2_transformed.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            s.defer(mincresample(img = anatomical_volume_2,
                                 xfm = xfm_concat,
                                 like = anatomical_volume_1_like,
                                 resampled = anatomical_volume_2_transformed,
                                 output_dir = output_dir))
            count_volume_2_transformed = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                    brain.name + "_" + count + "_stacked_2_transformed.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            s.defer(mincresample(img=count_volume_2,
                                 xfm=xfm_concat,
                                 like=count_volume_1_like,
                                 resampled=count_volume_2_transformed,
                                 output_dir=output_dir))

            #add the first section's like with the transformed second section
            s.defer(mincmath(imgs = [anatomical_volume_1_like, anatomical_volume_2_transformed],
                             result = anatomical_volume,
                             output_dir=output_dir))
            s.defer(mincmath(imgs=[count_volume_1_like, count_volume_2_transformed],
                             result=count_volume,
                             output_dir=output_dir))

#############################
# Step 4: Run autocrop to resample to isotropic
#############################
        anatomical_volume_isotropic = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                        brain.name + "_" + anatomical + "_stacked_isotropic.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
        anatomical_volume_isotropic_results = s.defer(autocrop(
            isostep = options.stacks_to_volume.plane_resolution,
            img = anatomical_volume,
            autocropped = anatomical_volume_isotropic
        ))

        count_volume_isotropic = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                 brain.name + "_" + count + "_stacked_isotropic.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
        count_volume_isotropic_results = s.defer(autocrop(
            isostep = options.stacks_to_volume.plane_resolution,
            img = count_volume,
            autocropped = count_volume_isotropic,
            nearest_neighbour = True
        ))

#############################
# Step 5: Run autocrop to pad the isotropic images
#############################
        x_pad = options.autocrop.x_pad
        y_pad = options.autocrop.y_pad
        z_pad = options.autocrop.z_pad

        anatomical_padded = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                              brain.name + "_" + anatomical + "_padded.mnc"))
        anatomical_pad_results = s.defer(autocrop(
            img = anatomical_volume_isotropic,
            autocropped = anatomical_padded,
            x_pad = x_pad,
            y_pad = y_pad,
            z_pad = z_pad
        ))
        all_anatomical_pad_results.append(anatomical_pad_results)
        reconstructed_mincs.append(anatomical_pad_results)

        count_padded = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                        brain.name + "_" + count + "_padded.mnc"))
        count_resampled = MincAtom(os.path.join(output_dir, pipeline_name + "_resampled",
                                          brain.name + "_" + count + "_resampled.mnc"))
        atlas_resampled = MincAtom(os.path.join(output_dir, pipeline_name + "_resampled",
                                          "atlas_to_" + brain.name + "_resampled.mnc"))
        count_pad_results = s.defer(autocrop(
            img = count_volume_isotropic,
            autocropped = count_padded,
            x_pad = x_pad,
            y_pad = y_pad,
            z_pad = z_pad
        ))
        all_count_pad_results.append(count_pad_results)
        all_count_resampled.append(count_resampled)
        all_atlas_resampled.append(atlas_resampled)
        reconstructed_mincs.append(count_pad_results)

    csv_file = pd.read_csv(options.application.csv_file)
    reconstructed = pd.DataFrame({'brain_directory': [brain.brain_directory.path for brain in brains],
                                  'z_slices': [brain.z for brain in brains],
                                  'z_resolution': [brain.z_resolution for brain in brains],
                                  'anatomical_padded': [anatomical_padded.path for anatomical_padded in all_anatomical_pad_results],
                                  'count_padded': [count_padded.path for count_padded in all_count_pad_results]})

    reconstructed = csv_file.merge(reconstructed)
    reconstructed.to_csv("reconstructed.csv", index=False)
    #TODO overlay them
    # s.defer(create_quality_control_images(imgs=reconstructed_mincs, montage_dir = output_dir,
    #     montage_output=os.path.join(output_dir, pipeline_name + "_stacked", "reconstructed_montage"),
    #                                       message="reconstructed_mincs"))

    s.defer(create_quality_control_images(imgs=all_anatomical_pad_results, montage_dir=output_dir,
                                          montage_output=os.path.join(output_dir, pipeline_name + "_stacked",
                                                                      "%s_montage" % anatomical),
                                          message="%s_mincs" % anatomical))
    s.defer(create_quality_control_images(imgs=all_count_pad_results, montage_dir=output_dir,
                                          montage_output=os.path.join(output_dir, pipeline_name + "_stacked",
                                                                      "%s_montage" % count),
                                          auto_range=True,
                                          message="%s_mincs" % count))
    return Result(stages=s, output=())

tv_recon_application = mk_application(parsers = [TV_stitch_parser,
                                                 cellprofiler_parser,
                                                 stacks_to_volume_parser,
                                                 autocrop_parser],
                                      pipeline = tv_recon_pipeline)

if __name__ == "__main__":
    tv_recon_application()
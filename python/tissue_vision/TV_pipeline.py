#!/usr/bin/env python3

import os
import sys

from configargparse import Namespace, ArgParser
from argparse import ArgumentError
from typing import Dict, List, Union

import pandas as pd

from pydpiper.core.stages import Stages, Result
from pydpiper.core.arguments import CompoundParser, AnnotatedParser, BaseParser
from pydpiper.core.util import maybe_deref_path
from pydpiper.execution.application import mk_application, execute
from pydpiper.core.files import FileAtom
from pydpiper.core.arguments import application_parser, registration_parser, execution_parser, parse
from pydpiper.minc.files import MincAtom
from pydpiper.minc.registration import autocrop, create_quality_control_images, check_MINC_input_files, lsq12_nlin, \
    get_linear_configuration_from_options, LinearTransType, get_nonlinear_component, concat_xfmhandlers, \
    get_registration_targets_from_init_model, xfmconcat
from pydpiper.pipelines.MBM import mbm, MBMConf, mk_mbm_parser
from pydpiper.pipelines.MAGeT import maget, maget_parsers, fixup_maget_options

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


def tissue_vision_pipeline(options):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name

    s = Stages()

    brains = get_brains(options.application) # List(Brain,...)

    # The new solution is to just wrap cellprofiler with cellprofiler.sh and load correct modules.
    env_vars = {}
    # env_vars['PATH'] = options.tissue_vision.cellprofiler.path
    # env_vars['PYTHONPATH'] = options.tissue_vision.cellprofiler.python_path
    # env_vars['LD_LIBRARY_PATH'] = options.tissue_vision.cellprofiler.ld_library_path
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
    reconstructed_mincs = []
    all_binary_resampled = []

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
                                                    TV_stitch_options = options.tissue_vision.TV_stitch,
                                                    Zstart=brain.z_start,
                                                    Zend=brain.z_end,
                                                    output_dir = output_dir
                                                    ))
            all_TV_stitch_results.append(TV_stitch_result)

        if brain.z_section:
            TV_stitch_result = s.defer(TV_stitch_wrap(brain_directory=brain.brain_directory,
                                                      brain_name=brain.name,
                                                      stitched=stitched[0 : brain.z_section - brain.z_start],
                                                      TV_stitch_options=options.tissue_vision.TV_stitch,
                                                      Zstart=brain.z_start,
                                                      Zend=brain.z_section - 1,
                                                      output_dir=output_dir
                                                      ))
            all_TV_stitch_results.append(TV_stitch_result)

            TV_stitch_result = s.defer(TV_stitch_wrap(brain_directory=brain.brain_directory,
                                                      brain_name=brain.name,
                                                      stitched=stitched[brain.z_section - brain.z_start:brain.z_end-1],
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
                                                  brain.name + "_" + anatomical + "_stacked.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
        binary_volume = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                      brain.name + "_" + binary + "_stacked.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))

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
                                                  brain.name + "_" + anatomical + "_stacked_1.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            smooth_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=smooths[0 : brain.z_section - brain.z_start],
                volume=smooth_volume_1,
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                uniform_sum=False,
                z_resolution=brain.z_resolution,
                output_dir=output_dir
            ))
            all_smooth_volume_results.append(smooth_slices_to_volume_results)

            smooth_volume_2 = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + anatomical + "_stacked_2.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            smooth_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=smooths[brain.z_section - brain.z_start:brain.z_end-1],
                volume=smooth_volume_2,
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                uniform_sum=False,
                z_resolution=brain.z_resolution,
                output_dir=output_dir
            ))
            all_smooth_volume_results.append(smooth_slices_to_volume_results)


            binary_volume_1 = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + binary + "_stacked_1.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            binary_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=binaries[0 : brain.z_section - brain.z_start],
                volume=binary_volume_1,
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                z_resolution=brain.z_resolution,
                uniform_sum=True,
                output_dir=output_dir
            ))
            all_binary_volume_results.append(binary_slices_to_volume_results)

            binary_volume_2 = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                  brain.name + "_" + binary + "_stacked_2.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            binary_slices_to_volume_results = s.defer(stacks_to_volume(
                slices=binaries[brain.z_section - brain.z_start:brain.z_end-1],
                volume=binary_volume_2,
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                z_resolution=brain.z_resolution,
                uniform_sum=True,
                output_dir=output_dir
            ))
            all_binary_volume_results.append(binary_slices_to_volume_results)

            #create minc slices for registration
            fixed_minc = s.defer(tif_to_minc(
                tif=smooths[brain.z_section - brain.z_start - 1],
                volume=MincAtom(smooths[brain.z_section - brain.z_start - 1].path.replace("tiff","mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked")),
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
                z_resolution=brain.z_resolution,
                output_dir=output_dir))
            moving_minc = s.defer(tif_to_minc(
                tif=smooths[brain.z_section - brain.z_start],
                volume=MincAtom(smooths[brain.z_section - brain.z_start].path.replace("tiff","mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked")),
                stacks_to_volume_options=options.tissue_vision.stacks_to_volume,
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
                                                         brain.name + "_" + anatomical + "_stacked_like.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            s.defer(get_like(img=smooth_volume_1, ref=smooth_volume_2,
                             like=smooth_volume_1_like, output_dir=output_dir))
            binary_volume_1_like = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                         brain.name + "_" + binary + "_stacked_like.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            s.defer(get_like(img=binary_volume_1, ref=binary_volume_2,
                             like=binary_volume_1_like, output_dir=output_dir))

            #transform the second section
            smooth_volume_2_transformed = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                brain.name + "_" + anatomical + "_stacked_2_transformed.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
            s.defer(mincresample(img = smooth_volume_2,
                                 xfm = xfm_concat,
                                 like = smooth_volume_1_like,
                                 resampled = smooth_volume_2_transformed,
                                 output_dir = output_dir))
            binary_volume_2_transformed = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                    brain.name + "_" + binary + "_stacked_2_transformed.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
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
                                                        brain.name + "_" + anatomical + "_stacked_isotropic.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
        smooth_volume_isotropic_results = s.defer(autocrop(
            isostep = options.tissue_vision.stacks_to_volume.plane_resolution,
            img = smooth_volume,
            autocropped = smooth_volume_isotropic
        ))
        all_smooth_volume_isotropic_results.append(smooth_volume_isotropic_results)

        binary_volume_isotropic = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                 brain.name + "_" + binary + "_stacked_isotropic.mnc"),
                                 output_sub_dir=os.path.join(output_dir, pipeline_name + "_stacked"))
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
        reconstructed_mincs.append(smooth_pad_results)

        binary_padded = MincAtom(os.path.join(output_dir, pipeline_name + "_stacked",
                                                        brain.name + "_" + binary + "_padded.mnc"))
        binary_resampled = MincAtom(os.path.join(output_dir, pipeline_name + "_resampled",
                                          brain.name + "_" + binary + "_resampled.mnc"))
        binary_pad_results = s.defer(autocrop(
            img = binary_volume_isotropic,
            autocropped = binary_padded,
            x_pad = x_pad,
            y_pad = y_pad,
            z_pad = z_pad
        ))
        all_binary_pad_results.append(binary_pad_results)
        all_binary_resampled.append(binary_resampled)
        reconstructed_mincs.append(binary_pad_results)

    csv_file = pd.read_csv(options.application.csv_file)
    reconstructed = pd.DataFrame({'brain_directory': [brain.brain_directory.path for brain in brains],
                                  'z_slices': [brain.z for brain in brains],
                                  'z_resolution': [brain.z_resolution for brain in brains],
                                  'smooth_padded': [smooth_padded.path for smooth_padded in all_smooth_pad_results],
                                  'binary_padded': [binary_padded.path for binary_padded in all_binary_pad_results]})

    reconstructed = csv_file.merge(reconstructed)
    reconstructed.to_csv("reconstructed.csv", index=False)
    #TODO overlay them
    # s.defer(create_quality_control_images(imgs=reconstructed_mincs, montage_dir = output_dir,
    #     montage_output=os.path.join(output_dir, pipeline_name + "_stacked", "reconstructed_montage"),
    #                                       message="reconstructed_mincs"))

    s.defer(create_quality_control_images(imgs=all_smooth_pad_results, montage_dir=output_dir,
                                          montage_output=os.path.join(output_dir, pipeline_name + "_stacked",
                                                                      "%s_montage" % anatomical),
                                          message="%s_mincs" % anatomical))
    s.defer(create_quality_control_images(imgs=all_binary_pad_results, montage_dir=output_dir,
                                          montage_output=os.path.join(output_dir, pipeline_name + "_stacked",
                                                                      "%s_montage" % binary),
                                          auto_range=True,
                                          message="%s_mincs" % binary))

#############################
# Step 6: Run MBM.py
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
# Step 7: Register consensus average to ABI tissuevision Atlas
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
# Step 8: Resample binary volumes to ABI tissuevision Atlas space
#############################
        all_full_xfms = []
        init_model = get_registration_targets_from_init_model(init_model_standard_file=options.mbm.lsq6.target_file,
                                                              output_dir=output_dir,
                                                              pipeline_name=pipeline_name)
        for mbm_xfm, binary_pad, binary_resampled in \
                zip(mbm_result.xfms.overall_xfm, all_binary_pad_results, all_binary_resampled):
            full_xfm = s.defer(xfmconcat([mbm_xfm.xfm, lsq12_nlin_result.xfm]))
            all_full_xfms.append(full_xfm)
            s.defer(mincresample(img=binary_pad,
                                 xfm=full_xfm,
                                 like=atlas_target,
                                 resampled=binary_resampled,
                                 output_dir=output_dir))

        transforms=transforms.assign(binary_resampled=[minc_atom.path for minc_atom in all_binary_resampled])
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
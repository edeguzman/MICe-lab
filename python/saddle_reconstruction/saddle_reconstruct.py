#!/usr/bin/env python3

""" Pipeline for reconstructing the data that comes from the MEMRI scan with the saddle coils using the varian scanner
The pipeline runs the following processes:
1) reconstruct mnc files from the fid files collected off the scanner using varian_recon.py
2) use lsq6 registration on the individual reps to create an average image for each coil
3) run distortion correction on the average images
4) crop the distortion corrected images

Created by: Elizabeth de Guzman
Date: April, 2016
"""

import os
import numpy

from configargparse import Namespace
from pydpiper.core.files import FileAtom
from pydpiper.core.arguments import lsq6_parser, CompoundParser, AnnotatedParser
from pydpiper.core.stages import Stages, Result
from pydpiper.execution.application import mk_application
from pydpiper.minc.files import MincAtom
from pydpiper.minc.registration import get_resolution_from_file, lsq6_nuc_inorm
from pydpiper.minc.registration import RegistrationTargets

from saddle_reconstruction.reconstruction import varian_recon_ge3dmice_saddle, dist_corr_saddle, crop_to_brain
from saddle_reconstruction.arguments      import varian_recon_parser, crop_to_brain_parser

def saddle_recon_pipeline(options):

    output_dir    = options.application.output_directory
    pipeline_name = options.application.pipeline_name
    fid_input_dir = options.saddle_recon.varian_recon.fid_input_directory

    # TODO check that the varian recon "output_file_name" doesn't name a directory, or if it does, that it matches the output directory

    #The input directory should contain a host of fid files that will be used for the reconstruction of the mnc files
    # TODO check that there are fid files in this directory
    # TODO check that all mandatory inputs are provided
    #Make list of input fid files, with location, and create a "FileAtom" type
    varian_fid_files = [fid_input_dir + "/fid" + str(num_fid)
                        for num_fid in range(0,options.saddle_recon.varian_recon.num_fids)]
    fids = [FileAtom(name) for name
            in varian_fid_files]

    # Varian recon will spit out images of the format <output_file_name>.<coil#>.<rep#>.mnc
    # TODO If .mnc isn't provided at the end of the output_file_name, then there is no "." before the coil number. Need to check and correct for this.
    # Al/ the files created will be spit out to the output_dir
    # Create list of files we expect to be produced by varian recon
    coil_list_0based = [int(x) for x in options.saddle_recon.varian_recon.mouse_list.split(',')]
    coil_list = [x+1 for x in coil_list_0based]
    file_name_base = os.path.splitext(options.saddle_recon.varian_recon.output_file_name)[0]
    varian_mnc_output = [os.path.join(output_dir, file_name_base + "." + str(coil) + "_" + str(rep) + ".mnc") for coil in coil_list for rep in range(0,options.saddle_recon.varian_recon.num_reps)]
    varian_coil_output = [str(coil) for coil in coil_list for rep in range(0,options.saddle_recon.varian_recon.num_reps)]
    recon_sub_dir = [file_name[:-6] + "_processed" for file_name in varian_mnc_output]
    imgs = [MincAtom(varian_mnc_output[k], pipeline_sub_dir=recon_sub_dir[k]) for k in range(0,len(varian_mnc_output))]

    s = Stages()

    #############################
    # Step 1: Run varian_recon.py
    #############################
    varian_recon_results = s.defer(varian_recon_ge3dmice_saddle(fids=fids,
                                                                imgs=imgs,
                                                                varian_recon_options=options.saddle_recon.varian_recon,
                                                                output_dir=output_dir))

    # Hold results obtained in the loop
    all_lsq6_results = []
    all_dc_results = []
    all_crop_results = []

    # Loop through all the coils
    for icoil in coil_list:
        icoil_imgs = numpy.array(imgs)[numpy.where(numpy.array(varian_coil_output) == str(icoil))[0]]
        icoil_varian_mnc_output = numpy.array(varian_mnc_output)[numpy.where(numpy.array(varian_coil_output) ==
                                                                             str(icoil))[0]]

        ###########################
        # Step 2: lsq6 registration
        ###########################
        # TODO change functionality of get_resolution_from_file so that it can be deferred
        lsq6_dir = os.path.join(output_dir, file_name_base + "." + str(icoil) + "_lsq6")
        target_dir = os.path.join(output_dir, file_name_base + "." + str(icoil) + "_target_file")

        resolution = options.registration.resolution
        #resolution = (options.registration.resolution or
        #              get_resolution_from_file(icoil_varian_mnc_output[0]))
        options.registration = options.registration.replace(resolution=resolution)

        target_file = MincAtom(name=icoil_varian_mnc_output[0], pipeline_sub_dir=target_dir)
        targets = RegistrationTargets(registration_standard=target_file,
                                   xfm_to_standard=None,
                                   registration_native=None)

        lsq6_result = s.defer(lsq6_nuc_inorm(imgs=icoil_imgs,
                                             resolution=resolution,
                                             registration_targets=targets,
                                             lsq6_dir=lsq6_dir,
                                             lsq6_options=options.saddle_recon.lsq6))
        all_lsq6_results.append(lsq6_result)

        ###########################
        # Step 3: distortion correct lsq6 output image
        ###########################
        lsq6_file = MincAtom(os.path.join(lsq6_dir,"average.mnc"), pipeline_sub_dir=lsq6_dir)
        dc_lsq6_file = MincAtom(os.path.join(lsq6_dir,"average.aug2015_dist_corr.mnc"), pipeline_sub_dir=lsq6_dir)
        dc_result = s.defer(dist_corr_saddle(img=lsq6_file,dc_img=dc_lsq6_file))
        all_dc_results.append(dc_result)

        ###########################
        # Step 4: crop distortion corrected lsq6 output image
        ###########################
        cropped_dc_lsq6_file = MincAtom(os.path.join(lsq6_dir,"average.aug2015_dist_corr.cropped.mnc"), pipeline_sub_dir=lsq6_dir)
        crop_result = s.defer(crop_to_brain(img=dc_lsq6_file,
                                            cropped_img=cropped_dc_lsq6_file,
                                            crop_to_brain_options=options.saddle_recon.crop_to_brain))
        all_crop_results.append(crop_result)

    return Result(stages=s, output=Namespace(varian_output=varian_recon_results,lsq6_output=all_lsq6_results,
                                             dc_output=all_dc_results,crop_output=all_crop_results))


saddle_recon_parser = CompoundParser([varian_recon_parser,lsq6_parser,crop_to_brain_parser])

saddle_recon_application = mk_application(parsers=[AnnotatedParser(parser=saddle_recon_parser, namespace='saddle_recon')],
                                          pipeline=saddle_recon_pipeline)

if __name__ == "__main__":
    saddle_recon_application()

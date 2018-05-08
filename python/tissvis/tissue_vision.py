#!/usr/bin/env python3

import os

from configargparse import Namespace
from typing import Dict, List

import pandas as pd

from pydpiper.core.stages import Stages, Result
from pydpiper.core.arguments import CompoundParser, AnnotatedParser
from pydpiper.execution.application import mk_application
from pydpiper.core.files import FileAtom, explode

from tissvis.arguments import TV_stitch_parser, cellprofiler_parser
from tissvis.reconstruction import TV_stitch_wrap, cellprofiler_wrap
from tissvis.TV_stitch import get_params

class Brain(object):
    def __init__(self,
                 brain_directory: FileAtom,
                 name: str,
                 slice_stitched: FileAtom = None,
                 slice_directory: FileAtom = None,
                 x: int = None,
                 y: int = None,
                 z: int = None) -> None:
        self.brain_directory = brain_directory
        self.name = name
        self.slice_stitched = slice_stitched
        self.slice_directory = slice_directory
        self.x = x
        self.y = y
        self.z = z


def get_brains(options):
    if options.files:
        raise ValueError("you used --files; please use --csv-file")

    csv = pd.read_csv(options.csv_file)

    brains = [Brain(FileAtom(brain_directory), brain_name)
              for brain_directory, brain_name in zip(csv.brain_directory, csv.brain_name)]
    return brains


def tissue_vision_pipeline(options):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name

    s = Stages()

    brains = get_brains(options.application) # List(Brain,...)

    # Hold results obtained in the loop
    all_TV_stitch_results = []
    all_cellprofiler_results = []

#############################
# Step 1: Run TV_stitch.py
#############################
    for brain in brains: #TODO remove this property
        brain.slice_directory = FileAtom(os.path.join(output_dir, pipeline_name + "_stitched", brain.name))

        stitched = []
        brain.x, brain.y, brain.z = get_params(os.path.join(brain.brain_directory.path, brain.name))
        for z in range (1, brain.z+1):
            brain.slice_stitched = FileAtom(os.path.join(brain.slice_directory.path, brain.name + "_Z%04d.tif" % z))
            stitched.append(brain.slice_stitched)

        #TODO remove this when you change the output to be proper
        if not (os.path.exists(brain.slice_directory.path)):
            os.makedirs(brain.slice_directory.path)

        TV_stitch_result = s.defer(TV_stitch_wrap(brain_directory = brain.brain_directory,
                                                brain_name = brain.name,
                                                stitched = stitched,
                                                application_options = options.application,
                                                TV_stitch_options = options.tissue_vision.TV_stitch,
                                                output_dir = output_dir
                                                ))
        all_TV_stitch_results.append(TV_stitch_result)

#############################
# Step 2: Run cellprofiler
#############################

    env_vars = {}
    env_vars['PYTHONPATH'] = options.tissue_vision.cellprofiler.python2_path
    cppline = FileAtom(options.tissue_vision.cellprofiler.cellprofiler_pipeline)

    for brain in brains: #TODO remove this property
        brain.cp_directory = os.path.join(output_dir, pipeline_name + "_cellprofiler", brain.name)
        brain.batch_data = FileAtom(os.path.join(brain.cp_directory,"Batch_data.h5"))

        overLays = []
        smooths = []
        microglias = []
        for z in range(1, brain.z + 1):
            brain.slice_overLay = FileAtom(
                os.path.join(brain.cp_directory, brain.name + "_Z%04d_overLay.tiff" % z))
            brain.slice_smooth = FileAtom(
                os.path.join(brain.cp_directory, brain.name + "_Z%04d_smooth.tiff" % z))
            brain.slice_microglia = FileAtom(
                os.path.join(brain.cp_directory, brain.name + "_Z%04d_microglia.tiff" % z))
            overLays.append(brain.slice_overLay)
            smooths.append(brain.slice_smooth)
            microglias.append(brain.slice_microglia)

        s.defer(cellprofiler_wrap(slice_directory = brain.slice_directory,
                                  cellprofiler_pipeline = cppline,
                                  batch_data = brain.batch_data,
                                  overLays = brain.slice_overLay,
                                  smooths = brain.slice_smooth,
                                  microglias = brain.slice_microglia,
                                  output_dir = output_dir,
                                  env_vars = env_vars
                                 ))


    return Result(stages=s, output=Namespace(TV_stitch_output=all_TV_stitch_results))



#############################
# Combine Parser & Make Application
#############################

tissue_vision_parser = CompoundParser([TV_stitch_parser, cellprofiler_parser])

tissue_vision_application = mk_application(
    parsers=[AnnotatedParser(parser=tissue_vision_parser, namespace='tissue_vision')],
    pipeline=tissue_vision_pipeline)

#############################
# Run
#############################

if __name__ == "__main__":
    tissue_vision_application()
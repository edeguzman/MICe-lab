#!/usr/bin/env python3

import os

from configargparse import Namespace
import pandas as pd

from pydpiper.core.stages import Stages, Result
from pydpiper.core.arguments import CompoundParser, AnnotatedParser
from pydpiper.execution.application import mk_application
from pydpiper.core.files import FileAtom


from tissvis.arguments import TV_stitch_parser
from tissvis.reconstruction import TV_stitch_wrap

class Brain(object):
    def __init__(self,
                 brain_directory: FileAtom,
                 name: str,
                 slice_directory: FileAtom = None) -> None:
        self.brain_directory = brain_directory
        self.name = name
        self.slice_directory = slice_directory


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


#############################
# Step 1: Run TV_stitch.py
#############################
    for brain in brains:

        brain.slice_directory = FileAtom(os.path.join(output_dir, pipeline_name + "_stitched", brain.name))

        if not (os.path.exists(brain.slice_directory.path)):
            os.makedirs(brain.slice_directory.path)

        TV_stitch_results = s.defer(TV_stitch_wrap(brain_directory = brain.brain_directory,
                                                   brain_name = brain.name,
                                                   slice_directory = brain.slice_directory,
                                                   application_options = options.application,
                                                   TV_stitch_options = options.tissue_vision.TV_stitch,
                                                   output_dir = output_dir
                                                   ))

    return Result(stages=s, output=Namespace(TV_stitch_output=TV_stitch_results))

#############################
# Combine Parser & Make Application
#############################

tissue_vision_parser = CompoundParser([TV_stitch_parser])

tissue_vision_application = mk_application(
    parsers=[AnnotatedParser(parser=tissue_vision_parser, namespace='tissue_vision')],
    pipeline=tissue_vision_pipeline)

#############################
# Run
#############################

if __name__ == "__main__":
    tissue_vision_application()
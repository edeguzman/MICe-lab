#!/usr/bin/env python3

import os

from configargparse import Namespace
from pydpiper.core.stages import Stages, Result
from pydpiper.core.arguments import lsq6_parser, CompoundParser, AnnotatedParser
from pydpiper.execution.application import mk_application

from tissvis.reconstruction import saddle
from tissvis.arguments import TV_stitch_parser

def tissue_vision_pipeline(options):

    output_dir    = options.application.output_directory
    pipeline_name = options.application.pipeline_name
    tif_input_dir = options.tissue_vision.TV_stitch.tif_input_directory

    s = Stages()

    #############################
    # Step 1: Run TV_stitch.py
    #############################
    TV_stitch_results = s.defer(saddle( fids=fids,
                                        imgs=imgs,
                                        TV_stitch_options=options.tissue_vision.TV_stitch,
                                        output_dir=output_dir))

    return Result(stages=s, output=Namespace(TV_stitch_output=TV_stitch_results,))


tissue_vision_parser = CompoundParser([TV_stitch_parser])

tissue_vision_application = mk_application(parsers=[AnnotatedParser(parser=tissue_vision_parser, namespace='tissvis')],
                                           pipeline=tissue_vision_pipeline)

if __name__ == "__main__":
    tissue_vision_application()

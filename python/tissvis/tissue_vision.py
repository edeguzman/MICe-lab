#!/usr/bin/env python3

import os

from configargparse import Namespace
from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.core.arguments import CompoundParser, AnnotatedParser
from pydpiper.execution.application import mk_application

from tissvis.arguments import TV_stitch_parser

def TV_stitch_cmd(TV_stitch_options, output_dir: str):
    stage = CmdStage(inputs=(), outputs=(),
                     cmd=['TV_stitch.py',
                          # TODO if app                       '--verbose' if TV_stitch_options.verbose else '',
                          '--skip_tile_match' if TV_stitch_options.skip_tile_match else '',
                          '--scaleoutput' if TV_stitch_options.scaleoutput else '',
                          '--Zstart', TV_stitch_options.Zstart,
                          '--Zend', TV_stitch_options.Zend,
                          os.path.join(TV_stitch_options.top_level_input_directory, TV_stitch_options.name),
                          os.path.join(output_dir, TV_stitch_options.name)])
    print(stage.render())
    stage.set_log_file(log_file_name=os.path.join(output_dir, "tissvis.log"))

    return Result(stages=Stages([stage]), output=())


def tissue_vision_pipeline(options):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name
    top_level_input_dir = options.tissue_vision.TV_stitch.top_level_input_directory
    slice_output_dir = options.tissue_vision.TV_stitch.slice_output_directory
    s = Stages()

    #############################
    # Step 1: Run TV_stitch.py
    #############################
    TV_stitch_results = s.defer(TV_stitch_cmd(TV_stitch_options=options.tissue_vision.TV_stitch, output_dir=output_dir))

    return Result(stages=s, output=Namespace(TV_stitch_output=TV_stitch_results, ))

tissue_vision_parser = CompoundParser([TV_stitch_parser])

tissue_vision_application = mk_application(
    parsers=[AnnotatedParser(parser=tissue_vision_parser, namespace='tissue_vision')],
    pipeline=tissue_vision_pipeline)

#############################
# Run
#############################

if __name__ == "__main__":
    tissue_vision_application()
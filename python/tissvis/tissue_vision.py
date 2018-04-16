#!/usr/bin/env python3

import os

from configargparse import Namespace
from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.core.arguments import CompoundParser, AnnotatedParser
from pydpiper.execution.application import mk_application

from tissvis.arguments import TV_stitch_parser

def TV_stitch_cmd(application_options, TV_stitch_options, output_dir: str):
    stage = CmdStage(inputs=(), outputs=(),
                     cmd=['TV_stitch.py', '--clobber', '--keeptmp',
                          '--verbose' if application_options.verbose else '',
                          '--skip_tile_match' if TV_stitch_options.skip_tile_match else '',
                          '--scaleoutput %s' % TV_stitch_options.scale_output if TV_stitch_options.scale_output else '',
                          '--Zstart %s' % TV_stitch_options.Zstart if TV_stitch_options.Zstart else '',
                          '--Zend %s' % TV_stitch_options.Zend if TV_stitch_options.Zend else '',
                          '--Ystart %s' % TV_stitch_options.Ystart if TV_stitch_options.Ystart else '',
                          '--Yend %s' % TV_stitch_options.Yend if TV_stitch_options.Yend else '',
                          '--Xstart %s' % TV_stitch_options.Xstart if TV_stitch_options.Xstart else '',
                          '--Xend %s' % TV_stitch_options.Xend if TV_stitch_options.Xend else '',
                          '--nogradimag' if TV_stitch_options.no_gradient_image else '',
                          '--use_positions_file %s' % TV_stitch_options.use_positions_file \
                              if TV_stitch_options.use_positions_file else '',
                          '--save_positions_file %s' % TV_stitch_options.save_positions_file \
                              if TV_stitch_options.save_positions_file else '',
                          '--overlapx %s' % TV_stitch_options.overlapx if TV_stitch_options.overlapx else '',
                          '--overlapy %s' % TV_stitch_options.overlapy if TV_stitch_options.overlapy else '',
                          '--channel %s' % TV_stitch_options.channel if TV_stitch_options.channel else '',
                          '--Zref %s' % TV_stitch_options.Zref if TV_stitch_options.Zref else '',
                          '--Zstack_pzIcorr' if TV_stitch_options.inormalize_piezo_stack else '',
                          '--fastpiezo' if TV_stitch_options.fast_piezo else '',
                          '--short' if TV_stitch_options.short_int else '',
                          #'--file_type %s' % TV_stitch_options.use_positions_file if TV_stitch_options.use_positions_file else '',
                          #'--TV_file_type %s' % TV_stitch_options.use_positions_file if TV_stitch_options.use_positions_file el
                          '--use_IM' if TV_stitch_options.use_imagemagick else '',
                          os.path.join(TV_stitch_options.top_level_input_directory, TV_stitch_options.brain),
                          os.path.join(output_dir, TV_stitch_options.brain)])
    print(stage.render())
    stage.set_log_file(log_file_name=os.path.join(output_dir, "tissvis.log"))

    return Result(stages=Stages([stage]), output=())

def tissue_vision_pipeline(options):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name
    top_level_input_dir = options.tissue_vision.TV_stitch.top_level_input_directory
    #TODO slice_output_dir = options.tissue_vision.TV_stitch.slice_output_directory
    s = Stages()

    #############################
    # Step 1: Run TV_stitch.py
    #############################
    TV_stitch_results = s.defer(TV_stitch_cmd(application_options=options.application, \
            TV_stitch_options=options.tissue_vision.TV_stitch, output_dir=output_dir))

    return Result(stages=s, output=Namespace(TV_stitch_output=TV_stitch_results, ))


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
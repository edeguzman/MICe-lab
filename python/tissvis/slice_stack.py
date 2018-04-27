#!/usr/bin/env python3

import os

from configargparse import ArgParser, Namespace
from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.execution.application import mk_application

from tissvis.arguments import slice_stack_parser

def slice_stack_cmd(application_options, slice_stack_options, output_dir: str):
    stage = CmdStage(inputs=(), outputs=(),
                     cmd=['stacks_to_volume.py',
                          '--input-resolution %s' % slice_stack_options.input_resolution,
                          '--output-resolution %s' % slice_stack_options.output_resolution,
                          '--slice-gap %s' % slice_stack_options.slice_gap,
                          '%s' % slice_stack_options.input_directory,
                          '%s' % slice_stack_options.output_directory])
    print(stage.render())
    #TODO since CmdStage.output==None, this line is needed for now...
    stage.set_log_file(log_file_name=os.path.join(output_dir, "stack.log"))

    return Result(stages=Stages([stage]), output=())

def slice_stack_pipeline(options):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name
    s = Stages()

    slice_stack_results = s.defer(slice_stack_cmd(application_options=options.application, \
            slice_stack_options=options.slice_stack, output_dir=output_dir))

    return Result(stages=s, output=Namespace(slice_stack_output=slice_stack_results, ))

cellprofiler_application = mk_application(parsers=[slice_stack_parser], pipeline=slice_stack_pipeline)

if __name__ == "__main__":
    cellprofiler_application()

# stacks_to_volume.py --input-resolution 0.00137 --output-resolution 0.025 --slice-gap 0.075
    # /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/cellprofiler_fl/*_smooth.tiff
    # /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/volume/cube.mnc
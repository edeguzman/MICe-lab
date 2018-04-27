#!/usr/bin/env python3

import os

from configargparse import ArgParser, Namespace
from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.execution.application import mk_application

from tissvis.arguments import slice_stack_parser

def slice_stack_cmd(application_options, slice_stack_options, output_dir: str, env_vars):
    stage = CmdStage(inputs=(), outputs=(),
                     cmd=['stacks_to_volume.py',
                          '--input-resolution %s' % slice_stack_options.input_resolution,
                          '--output-resolution %s' % slice_stack_options.output_resolution,
                          '--slice-gap %s' % slice_stack_options.slice_gap,
                          '%s' % slice_stack_options.output],
                     env_vars = env_vars)
    print(stage.render())
    #TODO since CmdStage.output==None, this line is needed for now...
    stage.set_log_file(log_file_name=os.path.join(output_dir, "stack.log"))

    return Result(stages=Stages([stage]), output=())

def slice_stack_pipeline(options):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name
    s = Stages()

    env_vars={}
    env_vars['PYTHONPATH']=options.cellprofiler.python2_path

    cellprofiler_results = s.defer(cellprofiler_batch(application_options=options.application, \
            cellprofiler_options=options.cellprofiler, env_vars = env_vars, output_dir=output_dir))

    return Result(stages=s, output=Namespace(cellprofiler_output=cellprofiler_results, ))

cellprofiler_application = mk_application(parsers=[cellprofiler_parser], pipeline=cellprofiler_pipeline)

if __name__ == "__main__":
    cellprofiler_application()

# stacks_to_volume.py --input-resolution 0.00137 --output-resolution 0.025 --slice-gap 0.075
    # /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/cellprofiler_fl/*_smooth.tiff
    # /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/volume/cube.mnc
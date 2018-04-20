#!/usr/bin/env python3

import os

from configargparse import ArgParser, Namespace
from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.core.arguments import AnnotatedParser, BaseParser, CompoundParser
from pydpiper.execution.application import mk_application

from tissvis.arguments import cellprofiler_parser

def cellprofiler_batch(application_options, cellprofiler_options, env_vars):
    stage = CmdStage(inputs=(), outputs=(),
                     cmd=['cellprofiler', '-c', '-r',
                          '-p %s' % cellprofiler_options.cellprofiler_pipeline,
                          '-i %s' % cellprofiler_options.image_directory,
                          '-o %s' % cellprofiler_options.output_directory],
                     env_vars = env_vars)
    print(stage.render())

    return Result(stages=Stages([stage]), output=())

def cellprofiler_pipeline(options):
    pipeline_name = options.application.pipeline_name
    s = Stages()

    #############################
    # Step 1: Run cellprofiler.py
    #############################
    env_vars={}
    env_vars['PYTHONPATH']=cellprofiler_options.python2_path

    cellprofiler_results = s.defer(cellprofiler_batch(application_options=options.application, \
            cellprofiler_options=options.cellprofiler.cellprofiler, env_vars = env_vars))
    TV_stitch_results = s.defer(TV_stitch_cmd(application_options=options.application, \
            TV_stitch_options=options.tissue_vision.TV_stitch, output_dir=output_dir, env_vars=env_vars))

    return Result(stages=s, output=Namespace(cellprofiler_output=cellprofiler_results, ))


#############################
# Make Parser & Application
#############################

cellprofiler_parser = CompoundParser([cellprofiler_parser])

cellprofiler_application = mk_application(
    parsers=[AnnotatedParser(parser=cellprofiler_parser, namespace='cellprofiler')],
    pipeline=cellprofiler_pipeline)

#cellprofiler_application = mk_application(parsers=[cellprofiler_parser], pipeline=cellprofiler_pipeline)
#lsq6_application = mk_application(parsers=[lsq6_parser], pipeline=lsq6_pipeline)
#############################
# Run
#############################

if __name__ == "__main__":
    cellprofiler_application()

# cellprofiler -c -r -p /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/FindMicroglia_Gauss.cppipe \
# -i /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/Salter_Microglia_GFP_SM1_28Feb13_2x2x2/ \
# -o /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/cellprofiler/

# cellprofiler -c -r -p /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/cellprofiler/Batch_data.h5 -f 1 -l 1

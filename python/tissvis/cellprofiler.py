#!/usr/bin/env python3

import os

from configargparse import ArgParser, Namespace
from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.core.arguments import AnnotatedParser, BaseParser
from pydpiper.execution.application import mk_application


def cellprofiler_batch(application_options, cellprofiler_options, env_vars = env_vars):
    stage = CmdStage(inputs=(), outputs=(),
                     cmd=['cellprofiler', '-c', '-r',
                          '-p %s' % cellprofiler_options.pipeline_filename,
                          '-i %s' % cellprofiler_options.image_directory,
                          '-o %s' % cellprofiler_options.skip_tile_match],
                     env_vars = env_vars)
    print(stage.render())

    return Result(stages=Stages([stage]), output=())

def cellprofiler_pipeline(options, cellprofiler_options):
    pipeline_name = options.application.pipeline_name
    s = Stages()

    #############################
    # Step 1: Run cellprofiler.py
    #############################
    env_vars={}
    env_vars['PYTHONPATH']=cellprofiler_options.python2_path

    cellprofiler_results = s.defer(cellprofiler_batch(application_options=options.application, \
            cellprofiler_options=options.cellprofiler, env_vars = env_vars))

    return Result(stages=s, output=Namespace(cellprofiler_output=cellprofiler_results, ))


#############################
# Make Parser & Application
#############################

def _mk_cellprofiler_parser():
    p = ArgParser(add_help=False)
    p.add_argument("-p", dest="pipeline_filename",
                   type=str,
                   default=None,
                   help="Load this pipeline file or project on startup")
    p.add_argument("-i", dest="image_directory",
                   type=str,
                   default=None,
                   help="Make this directory the default input folder")
    p.add_argument("-o", dest="image_directory",
                   type=str,
                   default=None,
                   help="Make this directory the default output folder")
    p.add_argument("-f", dest="first_image_set",
                   type=int,
                   default=None,
                   help="The one-based index of the first image set to process")
    p.add_argument("-l", dest="last_image_set",
                   type=int,
                   default=None,
                   help="The one-based index of the last image set to process")
    p.add_argument("--python2-path", dest="python2_path",
                   type=str,
                   default=None,
                   help="Cellprofiler is only compatible with python2. Set your python2 path using this flag.")
    return p

cellprofiler_parser = AnnotatedParser(parser=BaseParser(_mk_cellprofiler_parser(), "cellprofiler"),
                                   namespace="cellprofiler")

cellprofiler_application = mk_application(
    parsers=[AnnotatedParser(parser=cellprofiler_parser, namespace='cellprofiler')],
    pipeline=cellprofiler_pipeline)

#############################
# Run
#############################

if __name__ == "__main__":
    cellprofiler_application()

# cellprofiler -c -r -p /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/FindMicroglia_Gauss.cppipe \
# -i /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/Salter_Microglia_GFP_SM1_28Feb13_2x2x2/ \
# -o /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/cellprofiler/

# cellprofiler -c -r -p /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/cellprofiler/Batch_data.h5 -f 1 -l 1

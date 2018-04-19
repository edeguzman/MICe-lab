#!/usr/bin/env python3

import os

from configargparse import ArgParser, Namespace
from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.core.arguments import AnnotatedParser, BaseParser
from pydpiper.execution.application import mk_application


def cellprofiler_cmd(application_options, cellprofiler_options, output_dir: str):
    stage = CmdStage(inputs=(), outputs=(),
                     cmd=['cellprofiler', '-c', '-r',
                          '-p' if application_options.verbose else '',
                          '--skip_tile_match' if cellprofiler_options.skip_tile_match else '',
                          '--scaleoutput %s' % cellprofiler_options.scale_output if cellprofiler_options.scale_output else '',
                          os.path.join(cellprofiler_options.top_level_input_directory, cellprofiler_options.brain),
                          os.path.join(output_dir, cellprofiler_options.brain)])
    print(stage.render())
    stage.set_log_file(log_file_name=os.path.join(output_dir, "tissvis.log"))

    return Result(stages=Stages([stage]), output=())

def cellprofiler_pipeline(options):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name
    top_level_input_dir = options.tissue_vision.cellprofiler.top_level_input_directory
    s = Stages()

    #############################
    # Step 1: Run cellprofiler.py
    #############################
    cellprofiler_results = s.defer(cellprofiler_cmd(application_options=options.application, \
            cellprofiler_options=options.tissue_vision.cellprofiler, output_dir=output_dir))

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

# cellprofiler -p /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/FindMicroglia_Gauss.cppipe -i /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/Salter_Microglia_GFP_SM1_28Feb13_2x2x2/ -o /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/cellprofiler/ -c -r

# cellprofiler -p /hpf/largeprojects/MICe/nwang/Salter_Microglia_2x2x2/cellprofiler/Batch_data.h5 -c -r -f 1 -l 1

#!/usr/bin/env python3

import os

from configargparse import ArgParser, Namespace
from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.execution.application import mk_application

from tissvis.arguments import cellprofiler_parser

def cellprofiler_wrap(options):


    cellprofiler_results = s.defer(cellprofiler_batch(application_options=options.application,
            cellprofiler_options=options.cellprofiler, env_vars = env_vars, output_dir=output_dir))

    return Result(stages=s, output=Namespace(cellprofiler_output=cellprofiler_results, ))


#!/usr/bin/env python3

import os

from configargparse import Namespace, ArgParser
from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.core.arguments import BaseParser, CompoundParser, AnnotatedParser
from pydpiper.execution.application import mk_application
from pydpiper.core.util import NamedTuple


def TV_stitch_cmd(TV_stitch_options, output_dir: str):
    stage = CmdStage(inputs=(), outputs=(),
                     cmd=['TV_stitch.py',
                          # TODO if app                       '--verbose' if TV_stitch_options.verbose else '',
                          '--skip_tile_match' if TV_stitch_options.skip_tile_match else '',
                          '--scaleoutput', TV_stitch_options.scaleoutput,
                          '--Zstart', TV_stitch_options.Zstart,
                          '--Zend', TV_stitch_options.Zend,
                          os.path.join(TV_stitch_options.top_level_input_directory, TV_stitch_options.name),
                          os.path.join(output_dir, TV_stitch_options.name)])
    print(stage.render())
    stage.set_log_file(log_file_name=os.path.join(output_dir, "tissvis.log"))

    return Result(stages=Stages([stage]), output=())


def tissue_vision_pipeline(options: TVSTITCHConf):
    output_dir = options.application.output_directory
    pipeline_name = options.application.pipeline_name
    top_level_input_dir = options.tissue_vision.TV_stitch.top_level_input_dir
    slice_output_dir = options.tissue_vision.TV_stitch.slice_output_directory
    s = Stages()

    #############################
    # Step 1: Run TV_stitch.py
    #############################
    TV_stitch_results = s.defer(TV_stitch_cmd(TV_stitch_options=options.tissue_vision.TV_stitch, output_dir=output_dir))

    return Result(stages=s, output=Namespace(TV_stitch_output=TV_stitch_results, ))


#############################
# Parser
#############################
import pdb;

pdb.set_trace()


def _mk_TV_stitch_parser():
    p = ArgParser(add_help=False)
    p.add_argument("--skip-tile-match", dest="skip_tile_match",
                   type=bool,
                   action="store_true", default=False,
                   help="Skip tile matching and place tiles on perfect grid (for debugging)")
    p.add_argument("--scale-output", dest="scale_output",
                   type=int,
                   default=None,  # TODO raise an error when this isn't specified
                   help="Multiply slice images by this value before saving to file")
    p.add_argument("--Zstart", dest="Zstart",
                   type=int,
                   default=None,
                   help="Z start index")
    p.add_argument("--Zend", dest="Zend",
                   type=int,
                   default=None,
                   help="Z end index")
    p.add_argument("--top-level-input-directory", dest="top_level_input_directory",
                   type=str,
                   default='.',
                   help="Top level directory where folders containing tifs for each slice are contained")
    p.add_argument("--name", dest="name",
                   type=str,
                   default=None,
                   help="Name of the brain")
    return p


TVSTITCHConf = NamedTuple('TVSTITCHConf', [
    ('skip_tile_match', bool),
    ('scale_putput', int),
    ('Zstart', int),
    ('Zend', int),
    ('top_level_input_directory', str),
    ('name', str)
])


def to_TV_stitch_conf(TV_stitch_args: Namespace) -> TVSTITCHConf:
    return TVSTITCHConf(**TV_stitch_args.__dict__)


TV_stitch_parser = AnnotatedParser(parser=BaseParser(_mk_TV_stitch_parser(), "TV_stitch"),
                                   namespace="TV_stitch", cast=to_TV_stitch_conf)

tissue_vision_parser = CompoundParser([TV_stitch_parser])

tissue_vision_application = mk_application(
    parsers=[AnnotatedParser(parser=tissue_vision_parser, namespace='tissue_vision')],
    pipeline=tissue_vision_pipeline)

#############################
# Run
#############################

if __name__ == "__main__":
    tissue_vision_application()
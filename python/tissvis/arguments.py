from configargparse import ArgParser, Namespace
from pydpiper.core.util import NamedTuple
from pydpiper.core.arguments import BaseParser, AnnotatedParser

def _mk_TV_stitch_parser():
    p = ArgParser(add_help=False)
    p.add_argument("--hello-world", dest="hello_world",
                   type=bool,
                   default=False,
                   help="Prints Hello World")
    p.add_argument("--tif-input-directory", dest="tif_input_directory",
                   type=str,
                   default='.',
                   help="Directory where tif files are located")
    p.add_argument("--skip-tile-match", dest="skip_tile_match",
                   type=bool,
                   default=False,
                   help="Skip tile matching and place tiles on perfect grid (for debugging)")
    return p

TVSTITCHConf = NamedTuple('TVSTITCHConf', [('hello_world', bool),
                                                 ('tif_input_directory', str),
                                                 ('skip_tile_match', bool)]
                          )

def to_TV_stitch_conf(TV_stitch_args : Namespace) -> TVSTITCHConf:
    return TVSTITCHConf(**TV_stitch_args.__dict__)

TV_stitch_parser = AnnotatedParser(parser=BaseParser(_mk_TV_stitch_parser(), "TV_stitch"),
                                   namespace="TV_stitch", cast=to_TV_stitch_conf)

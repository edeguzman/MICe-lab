from configargparse import ArgParser, Namespace
from pydpiper.core.util import NamedTuple
from pydpiper.core.arguments import BaseParser, AnnotatedParser

def _mk_TV_stitch_parser():
    p = ArgParser(add_help=False)
    p.add_argument("--skip-tile-match", dest="skip_tile_match",
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
    ('scale_output', int),
    ('Zstart', int),
    ('Zend', int),
    ('top_level_input_directory', str),
    ('name', str)
])


def to_TV_stitch_conf(TV_stitch_args: Namespace) -> TVSTITCHConf:
    return TVSTITCHConf(**TV_stitch_args.__dict__)


TV_stitch_parser = AnnotatedParser(parser=BaseParser(_mk_TV_stitch_parser(), "TV_stitch"),
                                   namespace="TV_stitch", cast=to_TV_stitch_conf)
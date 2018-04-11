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
                   default=None,  # TODO raise a warning when this isn't specified
                   help="Multiply slice images by this value before saving to file")
    p.add_argument("--Zstart", dest="Zstart",
                   type=int,
                   default=None,
                   help="Z start index")
    p.add_argument("--Zend", dest="Zend",
                   type=int,
                   default=None,
                   help="Z end index")
    p.add_argument("--Ystart", dest="Ystart",
                   type=int,
                   default=None,
                   help="Y start index")
    p.add_argument("--Yend", dest="Yend",
                   type=int,
                   default=None,
                   help="Y end index")
    p.add_argument("--Xstart", dest="Xstart",
                   type=int,
                   default=None,
                   help="X start index")
    p.add_argument("--Xend", dest="Xend",
                   type=int,
                   default=None,
                   help="X end index")
    p.add_argument("--top-level-input-directory", dest="top_level_input_directory",
                   type=str,
                   default='.',
                   help="Top level directory where folders containing tifs for each slice are contained")
    p.add_argument("--brain", dest="brain",
                   type=str,
                   default=None,
                   help="Name of the brain")
    p.add_argument("--no-gradient-image", dest="no_gradient_image",
                   action="store_true", default=False,
                   help="Do not use gradient and raw image combined for correlation")
    p.add_argument("--save-positions-file", dest="save_positions_file",
                   type=str,
                   default=None,
                   help="Save the final positions to file (for subsequent use with --use-positions-file)")
    p.add_argument("--use-positions-file", dest="use_positions_file",
                   type=str,
                   default=None,
                   help="Use an existing positions file instead of generation positions from the input")
    p.add_argument("--overlapx", dest="overlapx",
                   type=int,
                   default=None,
                   help="% tile overlap in x direction")
    p.add_argument("--overlapy", dest="overlapy",
                   type=int,
                   default=None,
                   help="% tile overlap in y direction")
    p.add_argument("--channel", dest="channel",
                   type=str,
                   default=None,
                   help="channel to stitch")
    p.add_argument("--Zref", dest="Zref",
                   type=int,
                   default=None,
                   help="Z plane reference during tiling")
    p.add_argument("--inormalize-piezo-stack", dest="inormalize_piezo_stack",
                   action="store_true", default=False,
                   help="Intensity normalize piezo stacked images")
    p.add_argument("--fast-piezo", dest="fast_piezo",
                   action="store_true", default=False,
                   help="Piezo stack tiles are stored consecutively instead of plane-wise")
    p.add_argument("--short-int", dest="short_int",
                   action="store_true", default=False,
                   help="Write short int data to file instead of byte")
    p.add_argument("--use-imagemagick", dest="use_imagemagick",
                   action="store_true", default=False,
                   help="Use imagemagick for preprocessing (old behaviour)")

    return p

# TODO figure out with this with "cast=to_TV_stitch_conf" gives ValueError. Code can run without this though.
# TVSTITCHConf = NamedTuple('TVSTITCHConf', [
#     ('skip_tile_match', bool),
#     ('scale_output', int),
#     ('Zstart', int),
#     ('Zend', int),
#     ('top_level_input_directory', str),
#     ('slice_output_directory', str),
#     ('brain', str)
# ])
#
#
# def to_TV_stitch_conf(TV_stitch_args: Namespace) -> TVSTITCHConf:
#     return TVSTITCHConf(**TV_stitch_args.__dict__)


TV_stitch_parser = AnnotatedParser(parser=BaseParser(_mk_TV_stitch_parser(), "TV_stitch"),
                                   namespace="TV_stitch") #, cast=to_TV_stitch_conf)
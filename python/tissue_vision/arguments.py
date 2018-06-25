from configargparse import ArgParser, Namespace
from pydpiper.core.util import NamedTuple
from pydpiper.core.arguments import BaseParser, AnnotatedParser



#############################################

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
    p.add_argument("--brain", dest="brain",
                   type=str,
                   default=None,
                   help="Name of the brain")
    p.add_argument("--no-gradient-image", dest="no_gradient_image",
                   action="store_true", default=False,
                   help="Do not use gradient and raw image combined for correlation")
    p.add_argument("--keep-stitch-tmp", dest="keep_tmp",
                   action="store_true", default=False,
                   help="Keep temporary files from TV_stitch.")
    p.add_argument("--save-positions-file", dest="save_positions_file",
                   action="store_true", default=False,
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

TV_stitch_parser = AnnotatedParser(parser=BaseParser(_mk_TV_stitch_parser(), "TV_stitch"),
                                   namespace="TV_stitch")

#############################################

def _mk_cellprofiler_parser():
    p = ArgParser(add_help=False)
    p.add_argument("--cellprofiler-pipeline", dest="cellprofiler_pipeline",
                   type=str,
                   default=None,
                   help="Load this pipeline file or project on startup")
    p.add_argument("--image-directory", dest="image_directory",
                   type=str,
                   default=None,
                   help="Make this directory the default input folder")
    p.add_argument("--output-directory", dest="output_directory",
                   type=str,
                   default=None,
                   help="Make this directory the default output folder")
    p.add_argument("--first-image-set", dest="first_image_set",
                   type=int,
                   default=None,
                   help="The one-based index of the first image set to process")
    p.add_argument("--last-image-set", dest="last_image_set",
                   type=int,
                   default=None,
                   help="The one-based index of the last image set to process")
    p.add_argument("--cellprofiler-python-path", dest="python_path",
                   type=str,
                   default=None,
                   help="Set your PYTHONPATH environment variable for cellprofiler using this flag.")
    p.add_argument("--cellprofiler-path", dest="path",
                   type=str,
                   default=None,
                   help="Set your PATH environment variable for cellprofiler using this flag.")
    p.add_argument("--java-home", dest="java_home",
                   type=str,
                   default=None,
                   help="Set your JAVA_HOME environment variable for cellprofiler using this flag.")
    p.add_argument("--binary-name", dest="binary_name",
                   type=str,
                   default=None,
                   help="Specify the name of the binary images outputted by cellprofiler.")
    p.add_argument("--anatomical-name", dest="anatomical_name",
                   type=str,
                   default=None,
                   help="Specify the name of the anatomical images outputted by cellprofiler.")
    return p

cellprofiler_parser = AnnotatedParser(parser=BaseParser(_mk_cellprofiler_parser(), "cellprofiler"),
                                      namespace="cellprofiler")

#############################################

def _mk_stacks_to_volume_parser():
    p = ArgParser(add_help=False)
    p.add_argument("--input-resolution", dest="input_resolution",
                   type=float,
                   default=0.00137,
                   #TODO write help for standard and high resolution (get this from Dulcie)
                   help="")
    p.add_argument("--plane-resolution", dest="plane_resolution",
                   type=float,
                   default=None,
                   help="")
    return p

stacks_to_volume_parser = AnnotatedParser(parser=BaseParser(_mk_stacks_to_volume_parser(), "stacks_to_volume"),
                                      namespace="stacks_to_volume")

#############################################

def _mk_autocrop_parser():
    p = ArgParser(add_help=False)
    p.add_argument("--x-pad", dest="x_pad",
                   type=float,
                   default=0,
                   help="Single number in mm will be added to both sides")
    p.add_argument("--y-pad", dest="y_pad",
                   type=float,
                   default=0,
                   help="Single number in mm will be added to both sides")
    p.add_argument("--z-pad", dest="z_pad",
                   type=float,
                   default=0,
                   help="Single number in mm will be added to both sides")
    return p

autocrop_parser = AnnotatedParser(parser=BaseParser(_mk_autocrop_parser(), "autocrop"),
                                      namespace="autocrop")
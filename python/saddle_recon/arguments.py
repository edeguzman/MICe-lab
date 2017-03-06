from pkg_resources import get_distribution  # type: ignore
import copy
import os
import time
from configargparse import ArgParser, Namespace  # type: ignore
from typing import Any, Callable, List, Optional
from pydpiper.core.arguments import BaseParser, CompoundParser, AnnotatedParser
from pydpiper.core.util import AutoEnum, NamedTuple
from pydpiper.minc.registration import InputSpace, to_lsq6_conf, RegistrationConf, LSQ12Conf

def _mk_varian_recon_parser():
    p = ArgParser(add_help=False)
    p.add_argument("--fid-input-dir", dest="fid_input_directory",
                   type=str,
                   default='.',
                   help="Directory where fid files are located")
    p.add_argument("--petable", dest="petable",
                   type=str, default=None,
                   help="petable name and location if not specified by procpar file and current directory")
    p.add_argument("--mouse_list", dest="mouse_list",
                   type=str, default=None,
                   help="reconstruct list of mice instead of all mice (zero-based indexing, comma-delimited, no spaces)")
    p.add_argument("--fermi_ellipse", dest="fermi_ellipse",
                   action="store_true", default=False,
                   help="apply fermi ellipse filter")
    p.add_argument("--grappa_coil_groupings", dest="grappa_coil_groupings",
                   type=str, default=None,
                   help="groups of cross-talking coils (example: 0,2,5;1,3,4,6 )")
    p.add_argument("--procpar_file_name", dest="procpar_file_name",
                   type=str, default=None,
                   help="procpar file name (if not simply 'procpar')")
    p.add_argument("--output-file-name", dest="output_file_name",
                   type=str,
                   default='',
                   help="Name to give output files")
    p.add_argument("--num-reps", dest="num_reps",
                   type=int, default=None,
                   help="Tell code how many reps you expect varian_recon to spit out for each file")
    p.add_argument("--num-fids", dest="num_fids",
                   type=int, default=None,
                   help="Tell code how many fid files you are providing")
    p.add_argument("--final-dc-crop-output-dir", dest="final_dc_crop_output_dir",
                   type=str,
                   default='.',
                   help="Directory where you want the final dc and cropped images to be placed.")
    return p


VARIANRECONConf = NamedTuple('VARIANRECONConf', [('fid_input_directory', str),
                                                 ('petable', str),
                                                 ('mouse_list', str),
                                                 ('fermi_ellipse', bool),
                                                 ('grappa_coil_groupings', str),
                                                 ('procpar_file_name', str),
                                                 ('output_file_name', str),
                                                 ('num_fids', int),
                                                 ('num_reps', int),
                                                 ('final_dc_crop_output_dir', str)]

                             )

def to_varian_recon_conf(varian_recon_args : Namespace) -> VARIANRECONConf:
    return VARIANRECONConf(**varian_recon_args.__dict__)

varian_recon_parser = AnnotatedParser(parser=BaseParser(_mk_varian_recon_parser(), "varian_recon"), namespace="varian_recon", cast=to_varian_recon_conf)

def _mk_crop_to_brain_parser():
    p = ArgParser(add_help=False)
    p.add_argument("--crop_bbox_x", dest="bbox_x",
                   type=float, default=0,
                   help="length of bounding box in x direction (default units in pixels)")
    p.add_argument("--crop_bbox_y", dest="bbox_y",
                   type=float, default=0,
                   help="length of bounding box in y direction (default units in pixels)")
    p.add_argument("--crop_bbox_z", dest="bbox_z",
                   type=float, default=0,
                   help="length of bounding box in z direction (default units in pixels)")
    p.add_argument("--crop_buffer_z", dest="buffer_z",
                   type=float, default=0,
                   help="Add forced buffer in z direction (default units in pixels) (often the images sit too far forward)")
    p.add_argument("--crop_mm_units", action="store_true",
                   dest="mm_units", default=False,
                   help="Units of shift are in mm instead of pixels")
    return p

CROPTOBRAINConf = NamedTuple('CROPTOBRAINConf', [('bbox_x', float),
                                                 ('bbox_y', float),
                                                 ('bbox_z', float),
                                                 ('buffer_z', float),
                                                 ('mm_units', bool)]

                             )

def to_crop_to_brain_conf(crop_to_brain_args : Namespace) -> CROPTOBRAINConf:
    return CROPTOBRAINConf(**crop_to_brain_args.__dict__)

crop_to_brain_parser = AnnotatedParser(parser=BaseParser(_mk_crop_to_brain_parser(), "crop_to_brain"), namespace="crop_to_brain", cast=to_crop_to_brain_conf)

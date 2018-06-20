from pydpiper.core.files import FileAtom
from pydpiper.core.stages import Stages, CmdStage, Result
from pydpiper.core.util   import NamedTuple
from pydpiper.minc.files  import MincAtom, xfmToMinc
from pydpiper.minc.containers import XfmHandler
from pydpiper.minc.registration import concat_xfmhandlers, invert_xfmhandler

from typing import List, Optional, Tuple

import os


def varian_recon_ge3dmice_saddle(fids: List[FileAtom],
                                 imgs: List[MincAtom],
                                 varian_recon_options,
                                 output_dir: str):

    stage = CmdStage(inputs=tuple(fids), outputs=tuple(imgs), memory=25,
                     cmd=['mri_recon.py', 'ge3dmice_sg_cylA' , '--petable_ordered_pairs', '--grappa_coil_decouple',
                          '--outputreps', '--phasedriftcorr',
                          '--mouse_list',varian_recon_options.mouse_list,
                          '--petable',varian_recon_options.petable,
                          '--grappa_coil_groupings',varian_recon_options.grappa_coil_groupings,
                          '--procpar_file_name', varian_recon_options.procpar_file_name,
                          '--clobber', '--fermi_ellipse' if varian_recon_options.fermi_ellipse else '',
                          varian_recon_options.fid_input_directory,
                          os.path.join(output_dir,varian_recon_options.output_file_name)])
    print(stage.render())
    stage.set_log_file(log_file_name=os.path.join(output_dir,"varian_recon.log"))

    return Result(stages=Stages([stage]), output=imgs)

def dist_corr_saddle(img: MincAtom,
                     dc_img: MincAtom):

    stage = CmdStage(inputs=(img,), outputs=(dc_img,), memory=4,
                     cmd=['saddle_coil_distortion_correction_august_2015.pl',
                          '-spawn','-output-dir', img.pipeline_sub_dir,
                          img.path])
    print(stage.render())
    stage.set_log_file(log_file_name=os.path.join(img.pipeline_sub_dir,"dist_corr.log"))

    return Result(stages=Stages([stage]), output=dc_img)

def dist_corr_basket(img: MincAtom,
                     dc_img: MincAtom):

    stage = CmdStage(inputs=(img,), outputs=(dc_img,), memory=4,
                     cmd=['distortion_correction_september_2014.pl',
                          '-spawn','-output-dir', img.pipeline_sub_dir,
                          img.path])
    print(stage.render())
    stage.set_log_file(log_file_name=os.path.join(img.pipeline_sub_dir,"dist_corr.log"))

    return Result(stages=Stages([stage]), output=dc_img)

def crop_to_brain(img: MincAtom,
                  cropped_img: MincAtom,
                  crop_to_brain_options):

    stage = CmdStage(inputs=(img,), outputs=(cropped_img,), memory=4,
                     cmd=['crop_to_brain.py',
                          '--bbox_x',str(crop_to_brain_options.bbox_x),
                          '--bbox_y',str(crop_to_brain_options.bbox_y),
                          '--bbox_z',str(crop_to_brain_options.bbox_z),
                          '--intermediate_bbox=1.5',
                          '--buffer_z',str(crop_to_brain_options.buffer_z),
                          '--clobber',
                          img.path, cropped_img.path])
    print(stage.render())
    stage.set_log_file(log_file_name=os.path.join(img.pipeline_sub_dir,"crop_to_brain.log"))

    return Result(stages=Stages([stage]), output=cropped_img)

def shift_modify_header(img: MincAtom,
                        shifted_img: MincAtom,
                        newx: float,
                        newy: float,
                        newz: float):
    s = Stages()
    #Copy file to new location
    stage = CmdStage(inputs=(img,), outputs=(shifted_img,), memory=1,
                     cmd=['cp', img.path, shifted_img.path])
    print(stage.render())
    s.add(stage)
    #Alter header of copied image to shift
    xspace_start = 'xspace:start='+newx
    yspace_start = 'yspace:start='+newy
    zspace_start = 'zspace:start='+newz
    stage = CmdStage(inputs=(shifted_img,), outputs=(shifted_img,), memory=1,
                     cmd=['minc_modify_header','-dinsert',xspace_start,
                          '-dinsert',yspace_start,
                          '-dinsert',zspace_start,
                          shifted_img.path])
    print(stage.render())
    s.add(stage)
    #Alter header of copied image with header modification
    append_history_string = ':history= >>> copy and shift: '+ shifted_img.path +' to '+ shifted_img.path
    stage = CmdStage(inputs=(shifted_img,), outputs=(shifted_img,), memory=1,
                     cmd=['minc_modify_header','-sappend',
                          append_history_string,
                          shifted_img.path])
    print(stage.render())
    s.add(stage)
    return Result(stages=s, output=shifted_img)

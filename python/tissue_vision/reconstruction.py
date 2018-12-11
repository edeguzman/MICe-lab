import os

from typing import Dict, List

from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.core.files import FileAtom
from pydpiper.minc.files import MincAtom, XfmAtom
from pydpiper.core.util import NamedTuple

import pyminc.volumes.factory as pyminc

import numpy as np

def TV_stitch_wrap(brain_directory: FileAtom,
                   brain_name: str,
                   stitched: List[FileAtom],
                   TV_stitch_options,
                   Zstart: int,
                   Zend: int,
                   output_dir: str):
#TODO inputs should be tiles not just brain_directory
    stage = CmdStage(inputs=(brain_directory,), outputs=tuple(stitched),
                     cmd=['TV_stitch.py', '--clobber',
                          #'--verbose',
                          '--Zstart %s' % Zstart,
                          '--Zend %s' % Zend,
                          '--save_positions_file %s_positions.txt' %
                          os.path.join(stitched[0].dir, brain_name + '_Zstart' + str(Zstart))
                          if TV_stitch_options.save_positions_file else "",
                          '--keeptmp' if TV_stitch_options.keep_tmp else "",
                          '--scaleoutput %s' % TV_stitch_options.scale_output if TV_stitch_options.scale_output else '',
                          # '--skip_tile_match' if TV_stitch_options.skip_tile_match else '',
                          # '--Ystart %s' % TV_stitch_options.Ystart if TV_stitch_options.Ystart else '',
                          # '--Yend %s' % TV_stitch_options.Yend if TV_stitch_options.Yend else '',
                          # '--Xstart %s' % TV_stitch_options.Xstart if TV_stitch_options.Xstart else '',
                          # '--Xend %s' % TV_stitch_options.Xend if TV_stitch_options.Xend else '',
                          # '--nogradimag' if TV_stitch_options.no_gradient_image else '',
                          # '--use_positions_file %s' % TV_stitch_options.use_positions_file
                          # if TV_stitch_options.use_positions_file else '',
                          # '--overlapx %s' % TV_stitch_options.overlapx if TV_stitch_options.overlapx else '',
                          # '--overlapy %s' % TV_stitch_options.overlapy if TV_stitch_options.overlapy else '',
                          # '--channel %s' % TV_stitch_options.channel if TV_stitch_options.channel else '',
                          # '--Zref %s' % TV_stitch_options.Zref if TV_stitch_options.Zref else '',
                          # '--Zstack_pzIcorr' if TV_stitch_options.inormalize_piezo_stack else '',
                          # '--fastpiezo' if TV_stitch_options.fast_piezo else '',
                          # '--short' if TV_stitch_options.short_int else '',
                          # '--file_type %s' % TV_stitch_options.use_positions_file if TV_stitch_options.use_positions_file else '',
                          # '--TV_file_type %s' % TV_stitch_options.use_positions_file if TV_stitch_options.use_positions_file el
                          # '--use_IM' if TV_stitch_options.use_imagemagick else '',
                          os.path.join(brain_directory.path, brain_name),
                          os.path.join(stitched[0].dir, brain_name)],
                     log_file = os.path.join(output_dir, "TV_stitch.log"))

    return Result(stages=Stages([stage]), output=(stitched))

CellprofilerMemCfg = NamedTuple("CellprofilerMemCfg",
                            [('base_mem', float),
                             ('mem_per_size', float)])

default_cellprofiler_mem_cfg = CellprofilerMemCfg(base_mem=14, mem_per_size=2e-7)

def cellprofiler_wrap(stitched: List[FileAtom],
                       cellprofiler_pipeline: FileAtom,
                       batch_data: FileAtom,
                       overLays: List[FileAtom],
                       anatomicals: List[FileAtom],
                       count: List[FileAtom],
                       Zstart: int,
                       Zend: int,
                       output_dir: str):
    s = Stages()

    stage = CmdStage(inputs=(stitched+[cellprofiler_pipeline]), outputs=(batch_data,),
                     cmd=['cellprofiler.sh', '-c', '-r',
                          '-p %s' % cellprofiler_pipeline.path,
                          '-i %s' % stitched[0].dir,
                          '-o %s' % batch_data.dir],
                     log_file = os.path.join(output_dir,"cellprofiler.log"))
    s.add(stage)
    def set_memory(stage: CmdStage, mem_cfg: NamedTuple, z):
        img_size = os.stat(stitched[0].path).st_size
        stage.setMem(mem_cfg.base_mem + img_size * mem_cfg.mem_per_size)
    #cellprofiler's indexing starts at 1. so if Zstart=5, z=1 gives the 5th slice!
    for z in range (1, Zend + 2 - Zstart):
        stage = CmdStage(inputs=(batch_data,), outputs=(overLays[z-1], anatomicals[z-1], count[z-1]),
                         cmd=['cellprofiler.sh', '-c', '-r',
                              '-p %s' % batch_data.path,
                              '-f %s' % z,
                              '-l %s' % z],
                         log_file=os.path.join(output_dir, "cellprofiler.log"))


        #z=z is evaluated at this point, to avoid the problem of how python handles evironment.
        stage.when_runnable_hooks.append(lambda s, z=z:
                                         set_memory(s, default_cellprofiler_mem_cfg, z))

        s.add(stage)
    return Result(stages=s, output=(overLays, anatomicals, count))

def stacks_to_volume( slices: List[FileAtom],
                      volume: MincAtom,
                      z_resolution: float,
                      stacks_to_volume_options,
                      output_dir: str,
                      uniform_sum: bool = False):
    stage = CmdStage(inputs=tuple(slices), outputs=(volume,),
                     cmd=['stacks_to_volume.py',
                          '--input-resolution %s' % stacks_to_volume_options.input_resolution,
                          '--output-resolution %s' % stacks_to_volume_options.plane_resolution,
                          '--slice-gap %s' % z_resolution,
                          '--uniform-sum' if uniform_sum else '',
                          ' '.join(slice.path for slice in slices.__iter__()), #is this hacky or the right way?
                          '%s' % volume.path],
                     log_file=os.path.join(output_dir, "stacks_to_volume.log"))

    return Result(stages=Stages([stage]), output=(volume))

#refer to the link below for changing these parameters:
#https://github.com/ANTsX/ANTs/wiki/Anatomy-of-an-antsRegistration-call
def antsRegistration(fixed: MincAtom,
                     moving: MincAtom,
                     transform: XfmAtom,
                     output_dir: str,
                     warped: str = "Warped.nii.gz",
                     inversewarped: str = "InverseWarped.nii.gz",
                     dimensionality: int = 3):
#TODO warped and inversewarped output to the working directory
    stage = CmdStage(inputs=(fixed, moving), outputs=(transform,),
                     cmd = ['antsRegistration', '--verbose 1', '--float 0', '--minc',
                            '--dimensionality %s' % dimensionality,
                            '--output [%s,%s,%s]' % (transform.path.replace('0_GenericAffine.xfm', ''), warped, inversewarped),
                            '--interpolation Linear',
                            '--use-histogram-matching 0',
                            '--winsorize-image-intensities [0.01,0.99]',
                            '--initial-moving-transform [%s,%s,1]' % (fixed.path, moving.path), #1 indicates center of mass
                            '--transform Translation[0.1]',
                            '--metric MI[%s,%s,1,32,Regular,0.25]' % (fixed.path, moving.path),
                            '--convergence [1000x500x250x0,1e-6,10]',
                            '--shrink-factors 12x8x4x2',
                            '--smoothing-sigmas 4x3x2x1vox'],
                     log_file=os.path.join(output_dir, "join_sections.log")
                     )

    return Result(stages=Stages([stage]), output=(transform))

def antsApplyTransforms(img: FileAtom,
                        transform: XfmAtom,
                        transformed: FileAtom,
                        output_dir: str,
                        dimensionality: int = 2):

    stage = CmdStage(inputs=(img, transform), outputs=(transformed,),
                     cmd = ['antsApplyTransform', '--verbose',
                            '--dimensionality %s' % dimensionality,
                            '--input %s' % img,
                            '--reference-image %s' % img,
                            '--output %s' % transformed,
                            '--transform %s' % transform,
                            ],
                     log_file=os.path.join(output_dir, "join_sections.log"))
    return Result(stages=Stages([stage]), output=(transformed,))

def tif_to_minc(tif: FileAtom,
              volume: MincAtom,
              z_resolution: float,
              stacks_to_volume_options,
              output_dir: str,
              uniform_sum: bool = False):
    stage = CmdStage(inputs=(tif,), outputs=(volume,),
                     cmd=['stacks_to_volume.py',
                          '--input-resolution %s' % stacks_to_volume_options.input_resolution,
                          '--output-resolution %s' % stacks_to_volume_options.plane_resolution,
                          '--slice-gap %s' % z_resolution,
                          '--uniform-sum' if uniform_sum else '',
                          '%s' % tif.path,
                          '%s' % volume.path],
                     log_file=os.path.join(output_dir, "join_sections.log"))
    return Result(stages=Stages([stage]), output=(volume))

def get_like(img: MincAtom,
             ref: MincAtom,
             like: MincAtom,
             output_dir: str):

    stage = CmdStage(inputs=(img,ref), outputs=(like,),
                     cmd=['mincreshape', '-clobber',
                          '-start {start}',
                          '-count {count}',
                          img.path,
                          like.path],
                     log_file=os.path.join(output_dir, "join_sections.log"))

    def set_params(stage: CmdStage):
        img_pyminc = pyminc.volumeFromFile(img.path)
        ref_pyminc = pyminc.volumeFromFile(ref.path)
        count = np.maximum(ref_pyminc.data.count, img_pyminc.data.count)
        sum = ref_pyminc.data.count[0] + img_pyminc.data.count[0]
        for index, arg in enumerate(stage.cmd):
            stage.cmd[index] = arg.format(start = '0,0,0,', count = '%s,%s,%s' % (sum,count[1],count[2]))

    stage.when_runnable_hooks.append(lambda stage: set_params(stage))

    return Result(stages=Stages([stage]), output=(like))

def get_through_plane_xfm(img: MincAtom,
                          xfm: XfmAtom,
                          output_dir: str):
    stage = CmdStage(inputs=(img,), outputs=(xfm,),
                     cmd=['param2xfm', '-clobber', '-translation', '0', '{y}', '0',
                          xfm.path],
                     log_file=os.path.join(output_dir, "join_sections.log"))

    def set_params(stage:CmdStage):
        img_pyminc = pyminc.volumeFromFile(img.path)
        count = img_pyminc.data.count
        separations = img_pyminc.data.separations
        for index, arg in enumerate(stage.cmd):
            stage.cmd[index] = arg.format(y = count[0] * separations[0])

    stage.when_runnable_hooks.append(lambda stage: set_params(stage))

    return Result(stages=Stages([stage]), output=(xfm))

def concat_xfm(xfms: List[XfmAtom],
               outxfm: XfmAtom,
               output_dir: str):
    stage = CmdStage(inputs=tuple(xfms), outputs = (outxfm,),
                     cmd=['xfmconcat', '-clobber',
                          ' '.join(xfm.path for xfm in xfms.__iter__()),
                          outxfm.path],
                     log_file=os.path.join(output_dir, "join_sections.log"))
    return Result(stages=Stages([stage]), output=(outxfm))

# this was hacky
# def mincresample(img: MincAtom,
#                  xfm: XfmAtom,
#                  like: MincAtom,
#                  resampled: MincAtom,
#                  output_dir: str,
#                  invert_xfm: bool = False):
#     stage = CmdStage(inputs=(xfm, like, img), outputs=(resampled,),
#                      cmd=['mincresample', '-clobber',
#                           '-invert_transformation' if invert_xfm else '',
#                           '-transform %s' % xfm.path,
#                           '-like %s' % like.path, img.path, resampled.path],
#                      log_file = os.path.join(output_dir, "join_sections.log"))
#     return Result(stages=Stages([stage]), output=resampled)

def mincmath(imgs: List[MincAtom],
             result: MincAtom,
             output_dir: str):
    stage = CmdStage(inputs=tuple(imgs), outputs=(result,),
                     cmd=['mincmath', '-clobber', '-add',
                          ' '.join(img.path for img in imgs.__iter__()),
                          result.path],
                     log_file=os.path.join(output_dir, "join_sections.log"))
    return Result(stages=Stages([stage]), output=result)
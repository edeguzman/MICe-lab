import os

from typing import Dict, List

from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.core.files import FileAtom
from pydpiper.minc.files import MincAtom
from pyminc.volumes.factory import volumeFromFile
from pydpiper.core.util import NamedTuple


def TV_stitch_wrap(brain_directory: FileAtom,
                   brain_name: str,
                   stitched: List[FileAtom],
                   application_options,
                   TV_stitch_options,
                   output_dir: str):

    stage = CmdStage(inputs=(brain_directory,), outputs=tuple(stitched),
                     cmd=['TV_stitch.py', '--clobber', '--keeptmp',
                          '--verbose',
                          '--save_positions_file %s_positions.txt' % brain_name
                          if TV_stitch_options.save_positions_file else "",
                          '--scaleoutput %s' % TV_stitch_options.scale_output if TV_stitch_options.scale_output else '',
                          '--skip_tile_match' if TV_stitch_options.skip_tile_match else '',
                          '--Zstart %s' % TV_stitch_options.Zstart if TV_stitch_options.Zstart else '',
                          '--Zend %s' % TV_stitch_options.Zend if TV_stitch_options.Zend else '',
                          '--Ystart %s' % TV_stitch_options.Ystart if TV_stitch_options.Ystart else '',
                          '--Yend %s' % TV_stitch_options.Yend if TV_stitch_options.Yend else '',
                          '--Xstart %s' % TV_stitch_options.Xstart if TV_stitch_options.Xstart else '',
                          '--Xend %s' % TV_stitch_options.Xend if TV_stitch_options.Xend else '',
                          '--nogradimag' if TV_stitch_options.no_gradient_image else '',
                          '--use_positions_file %s' % TV_stitch_options.use_positions_file
                          if TV_stitch_options.use_positions_file else '',
                          '--overlapx %s' % TV_stitch_options.overlapx if TV_stitch_options.overlapx else '',
                          '--overlapy %s' % TV_stitch_options.overlapy if TV_stitch_options.overlapy else '',
                          '--channel %s' % TV_stitch_options.channel if TV_stitch_options.channel else '',
                          '--Zref %s' % TV_stitch_options.Zref if TV_stitch_options.Zref else '',
                          '--Zstack_pzIcorr' if TV_stitch_options.inormalize_piezo_stack else '',
                          '--fastpiezo' if TV_stitch_options.fast_piezo else '',
                          '--short' if TV_stitch_options.short_int else '',
                          # '--file_type %s' % TV_stitch_options.use_positions_file if TV_stitch_options.use_positions_file else '',
                          # '--TV_file_type %s' % TV_stitch_options.use_positions_file if TV_stitch_options.use_positions_file el
                          '--use_IM' if TV_stitch_options.use_imagemagick else '',
                          os.path.join(brain_directory.path, brain_name),
                          os.path.join(stitched[0].dir, brain_name)],
                     log_file = os.path.join(output_dir, "TV_stitch.log"))

    return Result(stages=Stages([stage]), output=(stitched))

CellprofilerMemCfg = NamedTuple("CellprofilerMemCfg",
                            [('base_mem', float),
                             ('mem_per_size', float)])

default_cellprofiler_mem_cfg = CellprofilerMemCfg(base_mem=1e-5, mem_per_size=1e-7)

def cellprofiler_wrap(stitched: List[FileAtom],
                       cellprofiler_pipeline: FileAtom,
                       batch_data: FileAtom,
                       overLays: List[FileAtom],
                       smooths: List[FileAtom],
                       binaries: List[FileAtom],
                       Zstart: int,
                       Zend: int,
                       output_dir: str,
                       env_vars: Dict[str, str]):
    s = Stages()

    stage = CmdStage(inputs=(stitched+[cellprofiler_pipeline]), outputs=(batch_data,),
                     cmd=['cellprofiler', '-c', '-r',
                          '-p %s' % cellprofiler_pipeline.path,
                          '-i %s' % stitched[0].dir,
                          '-o %s' % batch_data.dir],
                     log_file = os.path.join(output_dir,"cellprofiler.log"),
                     env_vars = env_vars)
    s.add(stage)

    def set_memory(stage: CmdStage, mem_cfg: NamedTuple, z):
        img_size = os.stat(stitched[z - 1].path).st_size
        stage.setMem(mem_cfg.base_mem + img_size * mem_cfg.mem_per_size)

    for z in range (Zstart, Zend + 1):
        stage = CmdStage(inputs=(batch_data,), outputs=(overLays[z-Zstart], smooths[z-Zstart], binaries[z-Zstart]),
                         cmd=['cellprofiler', '-c', '-r',
                              '-p %s' % batch_data.path,
                              '-f %s' % z,
                              '-l %s' % z],
                         log_file=os.path.join(output_dir, "cellprofiler.log"),
                         env_vars=env_vars)


        #z=z is evaluated at this point, to avoid the problem of how python handles evironment.
        stage.when_runnable_hooks.append(lambda s, z=z:
                                         set_memory(s, default_cellprofiler_mem_cfg, z))

        s.add(stage)
    return Result(stages=s, output=(overLays, smooths, binaries))

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

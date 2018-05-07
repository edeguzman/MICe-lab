import os

from pydpiper.core.stages import Stages, Result, CmdStage
from pydpiper.core.files import FileAtom


def TV_stitch_wrap(brain_directory: FileAtom,
                   brain_name: str,
                   slice_directory: FileAtom,
                   application_options,
                   TV_stitch_options,
                   output_dir: str):
    stage = CmdStage(inputs=(brain_directory,), outputs=(slice_directory,),
                     cmd=['TV_stitch.py', '--clobber', '--keeptmp',
                          '--verbose' if application_options.verbose else '',
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
                          os.path.join(slice_directory.path, brain_name)],
                     log_file=os.path.join(output_dir, "TV_stitch.log"))

    print(stage.render())

    return Result(stages=Stages([stage]), output=(slice_directory))

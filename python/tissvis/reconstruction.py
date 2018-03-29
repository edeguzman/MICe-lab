from pydpiper.core.stages import Stages, CmdStage, Result
import os

def saddle(tifs, imgs, TV_stitch_options, output_dir: str):

    stage = CmdStage(inputs=tuple(tifs), outputs=tuple(imgs),
                     cmd=['TV_stitch.py', '--verbose', '--skip_tile_match',
                          '--scaleoutput',TV_stitch_options.scaleoutput,
                          '--Zstart',TV_stitch_options.Zstart,
                          '--Zend', TV_stitch_options.Zend,
                          TV_stitch_options.tif_input_directory,
                          os.path.join(output_dir,TV_stitch_options.output_file_name)])
    print(stage.render())
    stage.set_log_file(log_file_name=os.path.join(output_dir,"tissvis.log"))

    return Result(stages=Stages([stage]), output=imgs)
#!/usr/bin/env python3
#

import argparse
import pandas as pd
import numpy as np
import cv2
import re
# from pathlib import Path
from fastai.vision import Path

# import argparse
#
# parser = argparse.ArgumentParser(description='Process some integers.')
# parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                     help='an integer for the accumulator')
# parser.add_argument('--sum', dest='accumulate', action='store_const',
#                     const=sum, default=max,
#                     help='sum the integers (default: find the max)')
#
# args = parser.parse_args()
# print(args.accumulate(args.integers))

parser = argparse.ArgumentParser(description='Perform a maximuim intensity projection on each stack'
                                             ' of piezo slices for the entire brain')
parser.add_argument("--input-dir", dest="input_dir", type=str,required=True)
parser.add_argument("--output-dir", dest="output_dir", type=str, required=True)
parser.add_argument("--name", dest="name", type=str, required=True)
args = parser.parse_args()

name = args.name
input_dir = Path(args.input_dir)
input_dir_ls = input_dir.ls()
output_dir = Path(args.output_dir)

slice_dirs = sorted([path for path in input_dir_ls if path.is_dir()])
tiles = [sorted([path for path in directory.ls() if ".tif" in str(path)],
               key=lambda tif: int(re.search(r"(?<=-)\d+(?=_)", str(tif)).group(0))) for directory in slice_dirs]
slice_mosaics = sorted([path for directory in slice_dirs for path in directory.ls() if
                 ".txt" in path.as_posix() and "Mosaic" in path.as_posix()])
mosaic = [path for path in input_dir_ls if
                 ".txt" in str(path) and "Mosaic" in str(path)][0]

values = pd.read_table(mosaic, sep = ":", names = ["parameter", "value"], error_bad_lines=False)
mrows = int(values[values.parameter=="mrows"].value.item())
mcolumns = int(values[values.parameter=="mcolumns"].value.item())
layers = int(values[values.parameter=="layers"].value.item())
sections = int(values[values.parameter=="sections"].value.item())

# TODO code for checking that things matchup as expected
# slice_dirs = [top_dir/(name + "-" +'{:04d}'.format(i))  for i in range(1,sections+1)]

output_slice_dirs = [output_dir/name/slice_dir.name for slice_dir in slice_dirs]
output_slice_mosaics = [output_slice_dir/slice_mosaic.name for
                       output_slice_dir, slice_mosaic in zip(output_slice_dirs, slice_mosaics)]

N = mrows * mcolumns
mip = [[output_slice_dirs[i]/tiles[i][j].name for j in range(N)] for i in range(sections)]

(pd.read_table(mosaic.as_posix(), header=None, squeeze=True)
     .str.replace("Zscan:1", "Zscan:0")
     .str.replace("layers:5", "layers:1")
     .to_csv(output_dir/name/mosaic.name, header=False, index=False))

for i in range(sections):
    print("Currently working on the tiles in " + output_slice_dirs[i].as_posix() + ".")
    output_slice_dirs[i].mkdir(parents=True, exist_ok=True)

# TODO dynamically reduce any number of piezo slices
    (pd.read_table(slice_mosaics[i].as_posix(), header=None, squeeze=True)
     .str.replace("Zscan:1", "Zscan:0")
     .str.replace("layers:5", "layers:1")
     .to_csv(output_slice_mosaics[i].as_posix(), header=False, index=False))


    for j in range(N):
        img = np.maximum.reduce([
            cv2.imread(tiles[i][j + 0 * N].as_posix(), cv2.COLOR_BGR2GRAY),
            cv2.imread(tiles[i][j + 1 * N].as_posix(), cv2.COLOR_BGR2GRAY),
            cv2.imread(tiles[i][j + 2 * N].as_posix(), cv2.COLOR_BGR2GRAY),
            cv2.imread(tiles[i][j + 3 * N].as_posix(), cv2.COLOR_BGR2GRAY),
            cv2.imread(tiles[i][j + 4 * N].as_posix(), cv2.COLOR_BGR2GRAY),
        ])
        cv2.imwrite(mip[i][j].as_posix(), img)
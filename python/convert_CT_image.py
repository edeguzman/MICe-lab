#!/usr/bin/env python3

import argparse
from argparse import Namespace
import os
import re
import subprocess
import sys


def parser():
  p = argparse.ArgumentParser("Convert a series of 2D .tif files from the Skyscan CT scanner to a MINC file,"
                              " given the log file from the scan and optionally a regex to filter tif files"
                              " in that directory")
  p.add_argument("--regex", type=lambda x: x.encode(),
                 help="regular expression (NOT shell glob) to match files to convert"
                      " (don't forget to escape special characters)")
  # p.add_argument("input_directory", help="directory containing files")
  p.add_argument("log_file", help="path to log file from CT scanner")
  p.add_argument("output_file", help="name of output file to create")
  p.add_argument("--print-files-only", action="store_true", default=False)
  p.add_argument("--print-files", action="store_true", default=False)
  return p


def same_dimensions(t1, t2):
  return t1.width == t2.width and t1.length == t2.length and t1.bit_depth == t2.bit_depth


class WrongDimensions(TypeError): pass


def tiff_info(tiff_file):
  info = subprocess.run(["tiffinfo", tiff_file], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout

  width  = re.search(rb"Image Width: (?P<width>\d+) ", info).group('width')
  length = re.search(rb"Image Length: (?P<length>\d+)", info).group('length')


  bit_depth = re.search(rb"Bits/Sample: (?P<bit_depth>\d+)", info).group('bit_depth')

  return Namespace(width=int(width), length=int(length), bit_depth=int(bit_depth))


def main(argv):
  args = parser().parse_args(argv[1:])

  d = os.path.dirname(args.log_file)
  all_files = subprocess.run(["ls", "-v", d], stdout=subprocess.PIPE).stdout.split(b"\n")

  # TODO match base of logfile.*\.tiff? instead of .tiff?
  tiff_files = [f for f in all_files if re.search(pattern=args.regex or b".tiff?", string=f)]

  if args.print_files:
    print("number of files: %d" % len(tiff_files))
    for f in tiff_files:
      print(f)

  if args.print_files_only:
    print("Printing files only ... done.")
    sys.exit()

  #tiff_infos = [tiff_info(t) for t in tiff_files]
  # hopefully the tiff conversion catches this itself:
  #if not all((same_dimensions(t, tiff_files[0]) for t in tiff_files));
  #  raise WrongDimensions("dimensions of files don't match")
  #dimensions = tiff_infos[0]

  dimensions = tiff_info(os.path.join(d.encode(), tiff_files[0]))

  print("dimensions: %s" % dimensions)

  with open(args.log_file, 'r') as log_file:
    contents = log_file.read()
    image_resolution = float(re.search("\nPixel Size.*=([\d.]+)", contents).group(1)) / 1000.0

  print("pixel size: %f" % image_resolution)

  if dimensions.bit_depth == 8:
    extra_args = ["-byte", "-unsigned", "-real_range", "0", "255"]
  elif dimensions.bit_depth == 16:
    extra_args = ["-short", "-unsigned", "-real_range", "0", "65535"]
  else:
    raise ValueError("don't know what to do with bit_depth %d" % dimensions.bit_depth)

  rawtominc_cmd = (["rawtominc", "-clobber", "-2", "-zxy",
                    "-xstep", str(image_resolution),
                    "-ystep", str(image_resolution),
                    "-zstep", str(image_resolution),
                    "-xstart", "0", "-ystart", "0", "-zstart", "0"]
                   + extra_args
                   + [args.output_file, str(len(tiff_files)), str(dimensions.length), str(dimensions.width)])

  # could use `p` as a context manager
  p = subprocess.Popen(args=rawtominc_cmd, stdin=subprocess.PIPE)
  for img in tiff_files:
    out = subprocess.run(["convert", "-quiet", os.path.join(d.encode(), img), "GRAY:-"],
                         stdout=subprocess.PIPE).stdout
    try:
      p.stdin.write(out)
    except BrokenPipeError:
      pass

  p.stdin.flush()
  p.stdin.close()
  

if __name__ == "__main__":
  main(sys.argv)

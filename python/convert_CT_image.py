#!/usr/bin/env python3

import argparse
from argparse import Namespace
import os
import re
import subprocess
import sys
import time


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
  p.add_argument("--check-file-transfer-complete", dest='check_transfer', action="store_true", default=False,
                 help="Indicate whether the .tif files you want to convert might still be in transfer. If "
                      "so, the program will first check that the number of .tif files is stable over time [not the default].")
  return p


def same_dimensions(t1, t2):
  return t1.width == t2.width and t1.length == t2.length and t1.bit_depth == t2.bit_depth


class WrongDimensions(TypeError): pass


def tiff_info(tiff_file):
  info = subprocess.run(["tiffinfo", tiff_file], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True).stdout

  width  = re.search(rb"Image Width: (?P<width>\d+) ", info).group('width')
  length = re.search(rb"Image Length: (?P<length>\d+)", info).group('length')


  bit_depth = re.search(rb"Bits/Sample: (?P<bit_depth>\d+)", info).group('bit_depth')

  return Namespace(width=int(width), length=int(length), bit_depth=int(bit_depth))


def main(argv):
  args = parser().parse_args(argv[1:])

  d = os.path.dirname(args.log_file)
  all_files = subprocess.run(["ls", "-v", d], stdout=subprocess.PIPE, check=True).stdout.split(b"\n")
  
  # The organization of files produced by the CT scanner is such that all TIFF slices pertaining to the scan are called:
  # {log_base}XXXXXXXX.tif where each X is an integer, we can use this information to find the correct tif files to work with
  log_no_extension = os.path.splitext(args.log_file)[0]
  log_base         = os.path.basename(log_no_extension)
  match_string     = log_base + "[0-9].*tiff?"
  tiff_files = [f for f in all_files if re.search(pattern=args.regex or bytes(match_string,'ascii'), string=f)]

  # We need to make sure that the transfer of TIFF files has completed (this program possibly is called in the middle
  # of a transfer, which would mean that only part of the file would be reconstructed). We count the number of TIFF files
  # related to this log file, and check again after 30 seconds until the TIFF count is stable.
  if args.check_transfer:
    prev_tif_count = -1
    cur_tif_count = len(tiff_files)
    print("Number of tiff files found: " + str(cur_tif_count))
    while prev_tif_count != cur_tif_count:
      # first wait a bit
      time.sleep(30)
      # read all files in the directory again:
      all_files = subprocess.run(["ls", "-v", d], stdout=subprocess.PIPE, check=True).stdout.split(b"\n")
      tiff_files = [f for f in all_files if re.search(pattern=args.regex or bytes(match_string,'ascii'), string=f)]
      prev_tif_count = cur_tif_count
      cur_tif_count = len(tiff_files)
      print("Number of tiff files found: " + str(cur_tif_count))

  if args.print_files or args.print_files_only:
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
                         stdout=subprocess.PIPE, check=True).stdout
    try:
      p.stdin.write(out)
    except BrokenPipeError:
      pass
  
  p.stdin.flush()
  p.stdin.close()
  p.wait()
  
  # add information about the scan (from the log file) to the header of the output file
  # first the entire content of the log file:
  log_file_handle  = open(args.log_file, 'r')
  log_file_lines = log_file_handle.readlines()
  minc_modify_header_cmd = ["minc_modify_header", "-sinsert", "CT:logfile=\"" + ''.join(log_file_lines) + "\"", args.output_file]
  subprocess.run(minc_modify_header_cmd)
  
  # also some specific parts:
  for line in log_file_lines:
    if "=" in line:
      # strip spaces/newline from both the beginning and end of the value after the equal sign:
      value = line.split("=")[1].strip()
      
      if "Camera binning" in line:
        subprocess.run(["minc_modify_header", "-sinsert", "acquisition:binning=" + value, args.output_file])
      elif "Image Rotation" in line:
        subprocess.run(["minc_modify_header", "-sinsert", "acquisition:imagerotation=" + value, args.output_file])
      elif "Exposure (ms)" in line:
        subprocess.run(["minc_modify_header", "-sinsert", "acquisition:exposure=" + value, args.output_file])
      elif "Rotation Step (deg)" in line:
        subprocess.run(["minc_modify_header", "-sinsert", "acquisition:rotationstep=" + value, args.output_file])
      elif "Use 360 Rotation" in line:
        subprocess.run(["minc_modify_header", "-sinsert", "acquisition:rotation360=" + value, args.output_file])
      elif "Source Voltage (kV)" in line:
        subprocess.run(["minc_modify_header", "-sinsert", "acquisition:sourcevoltage=" + value, args.output_file])
      elif "Source Current (uA)" in line:
        subprocess.run(["minc_modify_header", "-sinsert", "acquisition:sourcecurrent=" + value, args.output_file])
      elif "Vertical Object Position (mm)" in line:
        subprocess.run(["minc_modify_header", "-sinsert", "acquisition:verticalobjectposition=" + value, args.output_file])
      elif "Object to Source (mm)" in line:
        subprocess.run(["minc_modify_header", "-sinsert", "acquisition:objecttosource=" + value, args.output_file])
      elif "Camera to Source (mm)" in line:
        subprocess.run(["minc_modify_header", "-sinsert", "acquisition:cameratosource=" + value, args.output_file])
      elif "Filter=" in line:
        subprocess.run(["minc_modify_header", "-sinsert", "acquisition:filter=\"" + value + "\"", args.output_file])



if __name__ == "__main__":
  main(sys.argv)

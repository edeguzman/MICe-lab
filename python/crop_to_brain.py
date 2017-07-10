#!/usr/bin/env python3

from argparse import ArgumentParser
from time import time, ctime
import sys
import string, os
import numpy 
from numpy.fft import fft,fftshift,ifft
import pyminc.volumes.factory
import subprocess
import re
import multiprocessing
import types

object_dtype = "float"
minc_dtype = 'short'

program_name = 'crop_to_brain.py'
description = 'Takes an input image and iteratively crops the image smaller and smaller searching for the centre of mass until the final bounding box size has been achieved'

#----------------------------------------------------------------------
#-------Functions

def parse_file_name(path):
        file_name = path.split('/')[-1]
        file_folder = path.split(file_name)[0]
        file_end = file_name.split('.')[-1]
        file_beginning = file_name.split('.'+file_end)[0]
        return {'file_beginning':file_beginning, 'file_end':file_end, 'file_folder':file_folder}

#----------------------------------------------------------------------
#-------Top level program
if __name__ == '__main__':
        #Usage
        usage = "Usage: "+program_name+" [options] input.mnc output.mnc\n"+"   or  "+program_name+" --help"
        parser = ArgumentParser(usage=usage, description=description)
        parser.add_argument("input",help="input minc file (format: input.mnc)")
        parser.add_argument("output",help="output minc file (format: output.mnc)")
        parser.add_argument("--clobber", action="store_true", dest="clobber", default=0, help="overwrite output file")
        parser.add_argument("--bbox_x", dest="bbox_x", default=0, help="length of bounding box in x direction (default units in pixels)")
        parser.add_argument("--bbox_y", dest="bbox_y", default=0, help="length of bounding box in y direction (default units in pixels)")
        parser.add_argument("--bbox_z", dest="bbox_z", default=0, help="length of bounding box in z direction (default units in pixels)")
        parser.add_argument("--buffer_z", dest="buffer_z", default=0, help="Add forced buffer in z direction (default units in pixels) (often the images sit too far forward)")
        parser.add_argument("--intermediate_bbox", dest="bbox_int", default=1.5, help="size length of intermediate bbox with respect to final bbox (default:1.5)")
        parser.add_argument("--mm_units", action="store_true", dest="mm_units", default=False, help="Units of shift are in mm instead of pixels")
        args = parser.parse_args()
        #------------------------
        #Check inputs
        #Check that both an input and output file have been provided
        #Get input and output file
        input_file=args.input
        output_file=args.output
        print(input_file)
        print(output_file)
        input_file_pieces = parse_file_name(input_file)
        tmp_file = input_file_pieces['file_folder'] + input_file_pieces['file_beginning'] + '_croptobrain_tmp.' + input_file_pieces['file_end']
        #If output file exists check that the clobber option is set
        if not args.clobber and os.path.exists(output_file): raise SystemExit("The --clobber option is needed to overwrite an existing file.")
        if os.path.exists(tmp_file): raise SystemExit("Trying to make a temporary file. Why do you already have a file with the desired name? Groan.")
        #Create history stamp for printing to mnc file later 
        start_time = time()
        history = '>>> %s [%s]: %s\n' % (ctime(start_time), os.getcwd(), " ".join(sys.argv))
        print(history)
        #---------------------------
        #Start doing useful things
        #---------------------------
        #Get bounding box dimensions in pixels
        # Read in volume information
            # Doesn't put volume data to disk
        input_volume = pyminc.volumes.factory.volumeFromFile(input_file, dtype = object_dtype, readonly = True)
        if args.mm_units:
                volume_separation = input_volume.separations # mm/pixel in original mnc image
                bbox_x = int(float(args.bbox_x) / volume_separation[0])
                bbox_y = int(float(args.bbox_y) / volume_separation[1])
                bbox_z = int(float(args.bbox_z) / volume_separation[2])
                buffer_z = int(float(args.buffer_z) / volume_separation[2])
        else:
                bbox_x = int(float(args.bbox_x))
                bbox_y = int(float(args.bbox_y))
                bbox_z = int(float(args.bbox_z))
                buffer_z = int(float(args.buffer_z))

        bbox_int = float(args.bbox_int)
        
        #Large FOV
        print("Cropping to large FOV:...")
        bimodal_output = subprocess.check_output(['mincstats', '-biModalT', input_file], shell=False)
        bimodal_threshold = re.findall(r'\d+',str(bimodal_output))[0]
        com_output = subprocess.check_output(['mincstats',input_file ,'-CoM', '-floor', bimodal_threshold], shell=False)
        com = numpy.array(re.findall(r'\d+',str(com_output)))[[0,2,4]]
        start_string = str(int(int(com[0])-(bbox_int)*bbox_x/2)) + ',' +  str(int(int(com[1])-(bbox_int)*bbox_y/2)) + ',' + str(int(int(com[2])-(bbox_int)*bbox_z/2))
        count_string =  str(int((bbox_int)*bbox_x)) + ',' +   str(int((bbox_int)*bbox_y)) + ',' +  str(int((bbox_int)*bbox_z))
        print(['mincreshape', input_file, tmp_file, '-start', start_string, '-count', count_string])
        subprocess.call(['mincreshape', input_file, tmp_file, '-start', start_string, '-count', count_string], shell=False)
        #Small FOV
        print("Cropping to small FOV:...")
        bimodal_output = subprocess.check_output(['mincstats', '-biModalT', tmp_file], shell=False)
        bimodal_threshold = re.findall(r'\d+',str(bimodal_output))[0]
        com_output = subprocess.check_output(['mincstats',tmp_file ,'-CoM', '-floor', bimodal_threshold], shell=False)
        com = numpy.array(re.findall(r'\d+',str(com_output)))[[0,2,4]]
        start_string = str(int(int(com[0])-bbox_x/2)) + ',' +  str(int(int(com[1])-bbox_y/2)) + ',' + str(int(int(com[2])-bbox_z/2-25+buffer_z))
        count_string = str(bbox_x) + ',' +  str(bbox_y) + ',' + str(bbox_z)
        print(['mincreshape', tmp_file, output_file, '-start', start_string, '-count', count_string])
        subprocess.call(['mincreshape', tmp_file, output_file, '-start', start_string, '-count', count_string], shell=False)

        append_history_string = ':history= '+history
        subprocess.call(['minc_modify_header','-sappend',append_history_string,output_file], shell=False)
        subprocess.call(['rm',tmp_file],shell=False)




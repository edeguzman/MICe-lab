#!/usr/bin/env python3
#
# tifstackMIP.py
#
# A script for taking in group of tif images and calculating the MIP given a certain stack size
# Created July 24, 2017
# Elizabeth de Guzman
################################

from optparse import OptionParser, Option, OptionValueError
from time import time, ctime
from sys import argv
import string, os
from os import listdir
from os.path import isfile, join
import numpy
from PIL import Image

program_name = 'tifstackMIP.py'

description = """A script for taking in group of tif images and calculating the maximum imtensity projection (MIP) for each stack given a stack size provided by the user. The script assumes that the tif files are provided in numerical order of stack"""


class FatalError(Exception):
    def __init__(self,args=None):
        self.msg = args

#----------------------------------------------------------------------
# top level program

if __name__ == '__main__':

    usage = "Usage: "+program_name+" <input_directory> <output_file_basename>\n"+\
            "   or  "+program_name+" --help";
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--clobber", action="store_true", dest="clobber",
                       default=0, help="overwrite output file")
    parser.add_option("--stack_size", dest="stacksize", type=int,
                       default=None, help="Number of tif files in each stack")

    options, args = parser.parse_args()

    try:
        if len(args)<2: raise FatalError("Specify input directory and output file.")
        inputdirectory = args[-2]
        outputfile = args[-1]
        if not options.clobber and os.path.exists(outputfile):
            raise FatalError("The --clobber option is needed to overwrite an existing file.")
        if not options.stacksize or (options.stacksize == 1):
            raise SystemExit("Please provide a logical stack size.")
    except FatalError as e:
        print('Error(%s):' % program_name, e.msg)
        raise SystemExit

    #Find all files in the provided directory and sort them in descending order by name
    inputfiles = [f for f in os.listdir(inputdirectory) if os.path.isfile(os.path.join(inputdirectory, f))]
    sortedinputfiles = sorted(inputfiles)
    #Verify that the number of files is divisible by stacksize given
    try:
        if len(sortedinputfiles)%options.stacksize != 0: raise FatalError("Number of files provided not divisible by stack size")
    except FatalError as e:
        print('Error(%s):' % program_name, e.msg)
        raise SystemExit

    #Open the first file and determine the dimensions
    firstfile = Image.open(os.path.join(inputdirectory,sortedinputfiles[0]))
    fileshape = firstfile.size
    #Create empty array to fill with each stack
    stack_array = numpy.zeros((options.stacksize, fileshape[1], fileshape[0]), "float")
    #
    numMIPs = len(sortedinputfiles) // options.stacksize
    for thismip in range (0,numMIPs):
        for thisslice in range (0,options.stacksize):
            im = Image.open(os.path.join(inputdirectory,sortedinputfiles[thismip*options.stacksize+thisslice]))
            stack_array[thisslice,:,:] = numpy.array(im)
        stack_mip = stack_array.max(axis=0)
        image_mip = Image.fromarray(stack_mip.astype(numpy.uint8))
        image_mip.save(outputfile+"_"+str(thismip+1)+".tif")









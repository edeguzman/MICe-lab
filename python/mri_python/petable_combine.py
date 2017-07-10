#!/usr/bin/env python
#
# petable_combine.py
#
# a script to generate petables for varian
# acquisition and reconstruction
#
# Created May 24, 2013

from sys import argv
from sys import path as syspath
syspath.append('/home/bjnieman/source/vnmr')
from mri_python.varian_read_file import parse_petable_file
from mri_python.cyltable_gen import petable_file_output
import string
import os
import getopt
import re
from numpy import *
from numpy.random import shuffle
from optparse import OptionParser, Option, OptionValueError

program_name = 'petable_combine.py'

#________________________________________________________________________________________________




#________________________________________________________________________________________________

if __name__ == '__main__':

    usage = """%s <file1> ... <fileN> <output file>
   or  %s --help
   
%s is a script for combining table files into a single file.
"""
    usage = usage % ((program_name, )*3)

    parser = OptionParser(usage)
    parser.add_option("--clobber", action="store_true", dest="clobber",
                       default=0, help="overwrite output file")
    parser.add_option("--nf",type="int",dest="nf",
                       default=1, help="desired nf for each computed table")
    parser.add_option("--npetables",type="int",dest="npetables",
                       default=1, help="number of output tables")
    parser.add_option("--interleave", action="store_true", dest="interleave",
                       default=0, help="interleave tables rather than appending (must be tables of same length)")

    options, args = parser.parse_args()

    outputfile = args[-1]
    try:
        if not options.clobber and os.path.exists(outputfile):
            raise FatalError("The --clobber option is needed to overwrite an existing file.")
        pefile_list = args[0:-1]
    except FatalError as e:
        print('Error(%s):' % program_name, e.msg)
        raise SystemExit

    #i1=array((),int)
    #i2=array((),int)
    #for ctable in pefile_list:
    #    ci1 = parse_petable_file(ctable,"t1")
    #    ci2 = parse_petable_file(ctable,"t2")
    #    i1 = append(i1,ci1)
    #    i2 = append(i2,ci2)

    #read tables into list first
    i1list = []
    i2list = []
    for ctable in pefile_list:
        i1list.append(parse_petable_file(ctable,"t1"))
        i2list.append(parse_petable_file(ctable,"t2"))
    ntables = len(i1list)

    #combine by interleaving or by concatenating
    i1=array((),int)
    i2=array((),int)
    if (options.interleave): #interleave
        ntrstotal = len(i1list[0])*ntables//options.nf
        for j in range(ntrstotal):
            i1=append(i1,i1list[j%ntables][(j//ntables)*options.nf:(1+(j//ntables))*options.nf])
            i2=append(i2,i2list[j%ntables][(j//ntables)*options.nf:(1+(j//ntables))*options.nf])
    else: #concatenate
        for j in range(ntables):
            i1=append(i1,i1list[j])
            i2=append(i2,i2list[j])

    petable_file_output(i1,i2,options.nf,options.npetables,outputfile,linelim=options.nf)

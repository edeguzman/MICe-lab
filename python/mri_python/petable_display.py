#!/usr/bin/env python
#
# petable_display.py
#
# a script to display petables for varian
# acquisition and reconstruction
#

from sys import argv
from sys import path as syspath
syspath.append('/home/bjnieman/source/vnmr')
from mri_python.varian_read_file import parse_petable_file
import getopt
import re
from numpy import *
from mri_python.varian_read_file import parse_petable_file
from optparse import OptionParser, Option, OptionValueError
from pylab import plot,imshow,show,colorbar

program_name = 'petable_display.py'


#----------------------------------------------------------------------------
# define program specific exception
class FatalError(Exception):
    def __init__(self,args=None):
        self.msg = args


def display_petable(inputfile,etl=False,etl_lines=False,range_tuple=None):
    #read table file
    t1list=parse_petable_file(inputfile,'t1')
    t2list=parse_petable_file(inputfile,'t2')
    print('petable contains %d (t1) and %d (t2) elements'%(len(t1list),len(t2list)))
    if not (range_tuple is None):
        t1list=t1list[range_tuple[0]:range_tuple[1]]
        t2list=t2list[range_tuple[0]:range_tuple[1]] 
        print('...using elements %d to %d for display'%(range_tuple[0],range_tuple[1]))
    Nt1 = len(t1list)
    if (len(t2list)==0):
        nreps = Nt1//4
        t2list = [0 for j in t1list]
        for j in range(nreps-1):
            for k in range(Nt1):
                t2list.append(j)
                t1list.append(t1list[k])
    pemat=zeros((max(t2list)-min(t2list)+1,max(t1list)-min(t1list)+1),int)
    if not etl:
        pemat[array(t2list)-min(t2list),array(t1list)-min(t1list)]=arange(len(t1list))%Nt1
    else:
        pemat[array(t2list)-min(t2list),array(t1list)-min(t1list)]=1+arange(len(t1list))%etl
    imshow(pemat,interpolation='nearest'); colorbar()
    if ((etl_lines) and (etl)):
        #default: show 5 across the image
        for j in range(0,len(t1list)//etl,len(t1list)//(5*etl)):
            plot(t1list[j*etl+arange(etl)],t2list[j*etl+arange(etl)],'k-')
    show()
    return 0

def display_nacqs_petable(petable_name):
    enc_array_1 = parse_petable_file(petable_name,'t1')
    enc_array_2 = parse_petable_file(petable_name,'t2')
    enc_array_1 -= min(enc_array_1)
    enc_array_2 -= min(enc_array_2)
    nv = max(enc_array_1)+1
    nv2 = max(enc_array_2)+1
    nacq = zeros( (nv2,nv),int)
    for j in range(len(enc_array_1)):
        nacq[enc_array_2[j],enc_array_1[j]]+=1
    imshow(nacq,interpolation='nearest'); colorbar()
    show()
    return 0


#----------------------------------------------------------------------
# top level program

if __name__ == '__main__':

    usage = """%s <table file>
   or  %s --help
   
%s is a script for display petables.
"""
    usage = usage % ((program_name, )*3)

    parser = OptionParser(usage)
    parser.add_option("--nacqs", action="store_true", dest="nacqs",
                       default=0, help="display number of acquistions in repeated file instead of order")
    parser.add_option("--etl",type="int",dest="etl",
                       default=0, help="echo train length for echo order display")
    parser.add_option("--etl_lines", action="store_true", dest="etl_lines",
                       default=0, help="display some lines indicating trains when etl mode is used")
    parser.add_option("--disprange",type="string",dest="disprange",
                       default=None, help='range of indices to display (e.g.: "0,100")')
    
    options, args = parser.parse_args()

    inputfile = args[-1]

    try:
        if (len(args)<1):
                raise FatalError(usage)
    except FatalError as e:
        print('Error(%s):' % program_name, e.msg)
        raise SystemExit

    if (options.disprange is None):
        rangetuple=None
    else:
        r_words=options.disprange.split(',')
        rangetuple=(int(r_words[0]),int(r_words[1]))
    
    if (options.nacqs):
        display_nacqs_petable(inputfile)
    elif (options.etl):
        display_petable(inputfile,options.etl,options.etl_lines,range_tuple=rangetuple)
    else:
        display_petable(inputfile)


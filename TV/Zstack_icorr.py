
import os
import re
from numpy import *
from numpy.linalg import lstsq
from TV_stitch import TEMPDIRECTORY,program_name,run_subprocess

def IM_intensity_average(imgfile):
    cmdstr = "identify -verbose %s | grep mean "%imgfile
    cmdout = run_subprocess(cmdstr)
    Imean = float(cmdout.split()[1])
    return Imean

def intensity_normalize_Zstack(inputfilelist,outputfilelist):
    meanlist=[]
    for cfile in inputfilelist:
        meanlist.append( IM_intensity_average(cfile) )
    meanarray = array(meanlist)
    if (len(meanarray)>3): #HERE: fit to exponential ???
        A = ones( (len(meanarray),2) ); A[:,0]=arange(len(meanarray))
        x,resids,rank,s = lstsq(A,log(meanarray))
        scalearray = exp(-1*x[0]*arange(len(meanarray)))
    else:
        scalearray = meanarray[0]/meanarray
    for j in range(len(scalearray)):
        cmdstr = 'convert %s -evaluate multiply %f %s'%(inputfilelist[j],scalearray[j],outputfilelist[j])
        cmdout = run_subprocess(cmdstr)
    return 0


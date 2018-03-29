#!/usr/bin/python
#
# CGvesselfilter.py
#
# filter TV images for vessels
#
# Created March 2015

import sys
import os
import argparse
import exceptions
from numpy import *
from PIL import Image
from scipy.ndimage import gaussian_filter 
from scipy import interpolate


############################################################

class FatalError(exceptions.Exception):
    def __init__(self,args=None):
        self.msg = args


#maximum_filter at every pixel is too slow,
#do it at every size-th pixel and then interpolate
def spline_maxI_approx(inpimg,size=50):
    Nlim=(inpimg.shape[0]/size,inpimg.shape[1]/size)
    maxvals=zeros((Nlim[0],Nlim[1]),float)
    for j in range(Nlim[0]):
        for k in range(Nlim[1]):
            start0=j*size-size/2; start1=k*size-size/2
            end0=(j+1)*size-size/2+1; end1=(k+1)*size-size/2+1
            if (start0<0): start0=0
            if (start1<0): start1=1
            if (end0>inpimg.shape[0]): end0=inpimg.shape[0]
            if (end1>inpimg.shape[1]): end1=inpimg.shape[1]
            sortedI=sort(inpimg[start0:end0,start1:end1].flat)
            maxvals[j,k]=mean(sortedI[-int(0.05*size**2)::]) #average over top few instead of just taking max
    y=size*arange(maxvals.shape[0])
    x=size*arange(maxvals.shape[1])
    interpfunc=interpolate.interp2d(x,y,maxvals,kind='cubic')
    return interpfunc(arange(inpimg.shape[1]),arange(inpimg.shape[0]))
    

######################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="A script that filters TissueVision vessel images.")
    parser.add_argument("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    parser.add_argument("--filtersigma",type=float,dest="filtersigma",default=2.0,\
                      help="gaussian filter sigma")
    parser.add_argument("--high-pass-sigma",type=float,dest="high_pass_sigma",default=15,\
                      help="high pass filter gaussian sigma")
    parser.add_argument("--max-window",type=int,dest="max_window",default=30,\
                      help="window size for estimating maximum field correction")
    parser.add_argument("--max-floor",type=float,dest="max_floor",default=25,\
                      help="minimum value for max field estimate")
    parser.add_argument('inputfile')
    parser.add_argument('outputfile')

    args = parser.parse_args()

    try:
        if not args.clobber and os.path.exists(args.outputfile):
            raise FatalError, \
               "The --clobber option is needed to overwrite an existing file."
    except FatalError, e:
        print 'Error(%s):' % __file__, e.msg
        raise SystemExit

    
    #read image into array
    imgh=Image.open(args.inputfile)
    imgarr=asarray(imgh)

    #blur image
    blurimg = gaussian_filter(imgarr,args.filtersigma)

    #high pass filter image
    hpftemp = blurimg.astype(float)-gaussian_filter(blurimg.astype(float),args.high_pass_sigma)
    hpfimg = where(hpftemp<0,0,hpftemp).astype('uint8')

    #rescale hpfimg based on local max in image
    maxinterp=spline_maxI_approx(blurimg,size=args.max_window)
    outfullscale=hpfimg/where(maxinterp>args.max_floor,maxinterp,args.max_floor)
    outimg=where(outfullscale>0,where(outfullscale>1,1,outfullscale),0)

    #output the image
    outfh=Image.fromarray((255*outimg).astype('uint8'))
    outfh.save(args.outputfile)



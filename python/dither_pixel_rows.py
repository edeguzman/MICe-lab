#!/usr/bin/python
#
# dither_pixel_rows.py
#
# kluge to fix row-to-row mismatch in some TV tiles
#
# Created March, 2013

from sys import argv
import getopt
import exceptions
import os
from optparse import OptionParser, Option, OptionValueError
from PIL import Image
from scipy import *
from numpy import *
from scipy.signal import correlate2d

program_name = 'dither_pixel_rows.py'

#----------------------------------------------------------------------------
# define program specific exception
class FatalError(exceptions.Exception):
    def __init__(self,args=None):
        self.msg = args

#---------------------------------------------------------------------------

def pick_shifts(img,cwidth=25,rheight=10,kstart=18,kend=820):
    posshiftlist=[]
    for j in range(rheight,img.shape[0]-rheight-1,2): 
        c1 = img[j-rheight:j+rheight:2,:]
        c2 = img[j-rheight+1:j+rheight+1:2,:]
        for k in range(kstart,kend-cwidth,1):
            c1snip=c1[:,k:k+cwidth].copy()
            c2snip=c2[:,k-3:k+cwidth+3].copy()
            c1snip-=mean(c1snip); c1snip/=std(c1snip)
            c2snip-=mean(c2snip); c2snip/=std(c2snip)
            test=correlate2d(c2snip,c1snip,mode='valid')/(rheight*cwidth)
            bestind = argmax(test[0,:])
            if ((bestind>0) and (bestind<6) and (test[0,bestind]>0.2)):
                p=polyfit(arange(3)-1,test[0,bestind-1:bestind+2],2)
                pos_shift=bestind-3-0.5*p[1]/p[0]
                posshiftlist.append([j,k,pos_shift])
    posshifts=array(posshiftlist,float)
    return posshifts

from pylab import plot,show
def fit_shifts(posshifts,img,plotshifts=False,order=10):
    p = polyfit(posshifts[:,1],posshifts[:,2],order)
    x = arange(img.shape[1])
    xshifts = polyval(p,x)
    if (plotshifts):
        for j in range(0,832,10):
            flag=posshifts[:,0]==j
            plot(posshifts[flag,1],posshifts[flag,2],'.')
        plot(x,xshifts,'k-')
        show()
    return xshifts

def shift_img_pixels(img,xshifts):
    imgmod = img.copy()
    for j in range(0,img.shape[0]):
        if (j%2)==0:
            imgmod[j,:] = interp(arange(img.shape[1]),arange(img.shape[1])+0.5*xshifts,img[j,:])
        else:
            imgmod[j,:] = interp(arange(img.shape[1]),arange(img.shape[1])-0.5*xshifts,img[j,:])
    return imgmod


#---------------------------------------------------------------------------

if __name__ == '__main__':

    usage = """%s <inputfile> <outputfile>
   or  %s --help

%s is a kluge to fix row-to-row mismatch in some TV tiles based on correlation from 
   one row to the next.
"""
    usage = usage % ((program_name, )*3)

    parser = OptionParser(usage)
    parser.add_option("--order",type="int", dest="order",
                       default=10, metavar="<polynomial order>",help="Order of polynomial fit")
    parser.add_option("--plot_shifts_and_fit", action="store_true", dest="plot_shifts_and_fit",
                       default=0, help="Plot estimated shifts and fit result down the columns.")
    parser.add_option("--clobber", action="store_true", dest="clobber",
                       default=0, help="Overwrite existing output file if present.")
    options, args = parser.parse_args()

    try:
        if len(args)<2: raise FatalError, "No input file specified."
        inputfile = args[-2]
        outputfile = args[-1]
        if not os.path.exists(inputfile):
            raise FatalError, "Cannot find %s."%inputfile
        if (os.path.exists(outputfile) and not options.clobber):
            raise FatalError, "%s exists. Use --clobber to overwrite."%outputfile
    except FatalError, e:
        print 'Error(%s):' % program_name, e.msg
        raise SystemExit

    ifh = Image.open(inputfile)
    img = asarray( ifh.getdata() ,dtype=float)
    img.shape = ifh.size

    posshifts = pick_shifts(img)
    xshifts = fit_shifts(posshifts,img,plotshifts=options.plot_shifts_and_fit,order=options.order)
    imgmod = shift_img_pixels(img,xshifts)

    imgmodh = Image.fromarray(imgmod.astype(int16))
    imgmodh.save(outputfile)



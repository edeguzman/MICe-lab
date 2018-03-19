#!/usr/bin/python
#
# TV_stitch.py
#
# stitch XY planes of TV tiles
#
# Created February 2012

from sys import argv,path
import string
import os
import getopt
import exceptions
from optparse import OptionParser, Option, OptionValueError
import re
from PIL import Image
from numpy import *
from scipy.linalg import inv as inverse
import glob
from Zstack_icorr import *

program_name = 'TV_stitch.py'
TV_LORES = 1.37 #um
TV_HIRES = 0.55 #um

#----------------------------------------------------------------------------
# define program specific exception
class FatalError(exceptions.Exception):
    def __init__(self,args=None):
        self.msg = args

#---------------------------------------------------------------------------

def gen_tempfile(descrip_str,ftype_str):
    tempstr = "/tmp/" + program_name + '_' + descrip_str \
              + '_' + str(os.getpid()) + '.' + ftype_str
    return tempstr

def compose_img_file_list(inputdirectory,input_prefix,starts=[None,None,None],ends=[None,None,None],index_length=3,filetype="tif"):
    #determine starts and ends from the directory, use user input instead if requested
    pattern_str=inputdirectory.rstrip('/')+'/'+input_prefix
    for caxis in ['_Z','_Y','_X']:
        pattern_str+=caxis
        for j in range(index_length): pattern_str+="[0-9]"
    pattern_str+="."+filetype
    img_list=sorted(glob.glob(pattern_str))
    first_file=img_list[0]
    last_file=img_list[-1]
    restr=inputdirectory.rstrip('/')+'/'+input_prefix+'.*'
    if (starts[0]==None): starts[0]=int( re.search(restr+'_Z\d+',first_file).group(0)[-index_length:] )
    if (starts[1]==None): starts[1]=int( re.search(restr+'_Y\d+',first_file).group(0)[-index_length:] )
    if (starts[2]==None): starts[2]=int( re.search(restr+'_X\d+',first_file).group(0)[-index_length:] )
    if (ends[0]==None): ends[0]=int( re.search(restr+'_Z\d+',last_file).group(0)[-index_length:] )
    if (ends[1]==None): ends[1]=int( re.search(restr+'_Y\d+',last_file).group(0)[-index_length:] )
    if (ends[2]==None): ends[2]=int( re.search(restr+'_X\d+',last_file).group(0)[-index_length:] )
    print starts
    print ends
    #regenerate img_list in case user starts or ends specifies only a subset of images to work from
    img_list=[] 
    coord_list=[]
    for j in range(starts[0],ends[0]+1):
        for k in range(starts[1],ends[1]+1):
            for m in range(starts[2],ends[2]+1):
                cname = inputdirectory.rstrip('/')+'/'+input_prefix+\
                        '_Z'+str(j).zfill(index_length)+\
                        '_Y'+str(k).zfill(index_length)+\
                        '_X'+str(m).zfill(index_length)+\
                        '.'+filetype                     
                img_list.append(cname)
                coord_list.append([j,k,m])
    return img_list,array(coord_list,int)

def weighted_linear_least_squares(Amat,bcol,w):
    matrixmultiply=dot
    npts = len(Amat)
    W = identity(npts,dtype=float)
    for j in range(npts):
        W[j,j]=w[j]
    Bmat= bcol
    temp_val = matrixmultiply(W,Amat)
    temp_val = matrixmultiply(transpose(Amat),temp_val)
    temp_valB = matrixmultiply(W,bcol)
    temp_valB = matrixmultiply(transpose(Amat),temp_valB)
    x = matrixmultiply(inverse(temp_val),temp_valB)
    yresids = bcol - matrixmultiply(Amat,x)
    return x,yresids

def compute_offsets(img_list,coord_array,overlapx=20.0,overlapy=15.0,offsety=0.0,offsetx=0.0,Zref=-1,\
                    Cthresh=0.6,Dthresh=12,zsearch=20):
    #determine offsets of a grid of images with overlap at the edges based on CCimages (opencv correlative template matching)
    Nfiles=len(img_list)
    (matX,matY)=Image.open(img_list[0]).size
    offsetX=int(matX*(1.0-overlapx/100.0))
    offsetY=int(matY*(1.0-overlapy/100.0))
    positions = zeros((Nfiles,3),float)
    positions[:,0] = coord_array[:,0].astype(float)
    uniqueZ = unique(coord_array[:,0])
    Nz = uniqueZ.shape[-1]
    if (Zref<0): Zref=Nz/2                            #start in middle and work out by default (should really make this data driven)
    else: Zref-=1                                     #subtract 1 for 0 vs 1 based indexing
    sorted_uniqueZ = empty( (Nz,) , int)
    sorted_uniqueZ[0:Nz-Zref] = uniqueZ[-(Nz-Zref):]
    sorted_uniqueZ[(Nz-Zref):] = uniqueZ[0:Zref][::-1]  
    searchx_dict={'x': 0.8*matX*overlapx/100.0,      #this assumes overlap is on the small side (~20% and not near 100%)
                  'y': 0.8*matY*overlapy/100.0,
                  'z': zsearch}
    searchy_dict={'x': 0.8*matX*overlapx/100.0,
                  'y': 0.8*matY*overlapy/100.0,
                  'z': zsearch}
    templx_dict={'x': 0.5*matX*overlapx/100.0, 'y': matY-2*searchy_dict['x'], 'z': matX-2*zsearch}
    temply_dict={'x': matX-2*searchx_dict['y'], 'y': 0.5*matY*overlapy/100.0, 'z': matY-2*zsearch}
    grid_offset_results=array([[0,0,0,0]],float)
    refPmatch_results=[]
    for zindex,currz in enumerate(sorted_uniqueZ):
        cz_inds = nonzero(coord_array[:,0]==currz)[0]
        XY_img_list = [img_list[k] for k in cz_inds]
        cz_nfiles = len(XY_img_list)
        Alist=[]    #coefficients of positions to produce offsets (i.e. 0 0 0 ... 1 -1 0 ... 0)
        blist=[]    #offset results
        wlist=[]    #stores CCimages correlation results, after **2 serves as weights for least squares fit
        zposlist=[] #the x or y position from the previous z-slice (NOT the z-position)
        adjustedoffsets=zeros( (cz_nfiles,2*coord_array.shape[1]) , float)
        for caxis in ['x','y','z']:
            for j in range(0,cz_nfiles):
                cfile=XY_img_list[j]
                #identify current reference image
                relpos = coord_array[cz_inds[j],:]-coord_array[:,:]
                if (caxis=='x'):
                    testpos = (relpos[:,0]==0)*(relpos[:,1]==0)*(relpos[:,2]==1)
                elif (caxis=='y'):
                    testpos = (relpos[:,0]==0)*(relpos[:,1]==1)*(relpos[:,2]==0)
                elif (caxis=='z'):
                    if ( (currz-sorted_uniqueZ[0:zindex])==1 ).any():
                        testpos = (relpos[:,0]==1)*(relpos[:,1]==0)*(relpos[:,2]==0)
                    elif ( (currz-sorted_uniqueZ[0:zindex])==-1 ).any():
                        testpos = (relpos[:,0]==-1)*(relpos[:,1]==0)*(relpos[:,2]==0)
                    else:
                        testpos = zeros((relpos.shape[0],),int)
                if not (testpos.any()):
                    continue
                refind = nonzero(testpos)[0][0] 
                reffile=img_list[refind]
                #run CCimages
                searchx=int(searchx_dict[caxis])
                searchy=int(searchy_dict[caxis])
                templx =int(templx_dict[caxis])
                temply =int(temply_dict[caxis])
                xoff   =(coord_array[cz_inds[j],2]-coord_array[refind,2])*offsetX+\
                        (coord_array[cz_inds[j],1]-coord_array[refind,1])*offsetx
                yoff   =(coord_array[cz_inds[j],1]-coord_array[refind,1])*offsetY+\
                        (coord_array[cz_inds[j],2]-coord_array[refind,2])*offsety
                newoffsets,CCresult,Imean=run_CCimages([reffile,cfile],search_width=searchx,\
                                                                       search_height=searchy,\
                                                                       xoff=xoff,\
                                                                       yoff=yoff,\
                                                                       templ_width=templx,\
                                                                       templ_height=temply)
                print '\n'
                #put results into Alist, blist, wlist and zposlist
                Arow=zeros((2*cz_nfiles,),float)
                Arow[2*j]=1.0 
                if not (caxis=='z'):
                    Arow[2*where(cz_inds==refind)[0][0]]=-1.0
                else:
                    newoffsets[1]+=positions[refind,1]
                    zposlist.append(positions[refind,1])
                Alist.append(Arow); blist.append(newoffsets[1]); wlist.append(CCresult) #CCresult*sqrt(Imean)
                Arow=zeros((2*cz_nfiles,),float)
                Arow[2*j+1]=1.0
                if not (caxis=='z'):
                    Arow[2*where(cz_inds==refind)[0][0]+1]=-1.0
                else:
                    newoffsets[0]+=positions[refind,2]
                    zposlist.append(positions[refind,2])
                Alist.append(Arow); blist.append(newoffsets[0]); wlist.append(CCresult) #CCresult*sqrt(Imean)
        #generate arrays for lsq
        Aarray = array(Alist,float)
        if (zindex==0): #for first slice, need to define a reference arbitrarily, othinds indexes into cz_inds
            refind = nonzero( (coord_array[cz_inds,1]==min(coord_array[cz_inds,1]))*\
                              (coord_array[cz_inds,2]==min(coord_array[cz_inds,2])) )[0]
            othinds = nonzero( (coord_array[cz_inds,1]!=min(coord_array[cz_inds,1]))+\
                               (coord_array[cz_inds,2]!=min(coord_array[cz_inds,2])) )[0]
            Aarray = delete(Aarray,refind,axis=1) #these lines assign position (0,0) to image at Z(Nz/2)_Y001_X001-->need a full cartesian grid
            Aarray = delete(Aarray,refind,axis=1)        
        else:
            othinds = arange(len(cz_inds))
        Barray = array(blist,float)
        warray = array(wlist,float)
        #if individual offsets are way off or correlation is really bad, replace offset result with guess based on the rest of the grid
        if (cz_nfiles>8):
            refPmatch_results.append(median(warray)) #refPmatch=sort(warray)[3*cz_nfiles/2]
            Pthresh=Cthresh*median(refPmatch_results)             #Pthresh=Cthresh*refPmatch
            #print "Discard threshold set to %f"%Pthresh
            #print "%d matches below threshold."%sum(less(warray[::2],Pthresh))
            grid_x = max(coord_array[cz_inds,2])-min(coord_array[cz_inds,2])+1
            grid_y = max(coord_array[cz_inds,1])-min(coord_array[cz_inds,1])+1
            yyinds = arange(2*(cz_nfiles-grid_y),2*(2*cz_nfiles-grid_x-grid_y),2)
            adjyy = median(take(Barray[yyinds],nonzero( greater(warray[yyinds],Pthresh) )[0]))
            yxinds = arange(1+2*(cz_nfiles-grid_y),2*(2*cz_nfiles-grid_x-grid_y),2) 
            adjyx = median(take(Barray[yxinds],nonzero( greater(warray[yxinds],Pthresh) )[0]))
            xyinds = arange(0,2*(cz_nfiles-grid_y),2)
            adjxy = median(take(Barray[xyinds],nonzero( greater(warray[xyinds],Pthresh) )[0]))
            xxinds = arange(1,2*(cz_nfiles-grid_y),2)
            adjxx = median(take(Barray[xxinds],nonzero( greater(warray[xxinds],Pthresh) )[0]))
            grid_offset_results = append(grid_offset_results,array([[adjxy,adjxx,adjyy,adjyx]],float),axis=0)
            Bideal = zeros(Barray.shape,float)
            Bideal[xyinds]=median(grid_offset_results[1:,0]) #adjxy
            Bideal[xxinds]=median(grid_offset_results[1:,1]) #adjxx
            Bideal[yyinds]=median(grid_offset_results[1:,2]) #adjyy
            Bideal[yxinds]=median(grid_offset_results[1:,3]) #adjyx
            if (len(zposlist)>0):
                Bideal[-len(zposlist):]=array(zposlist,float)
            Barray[:]=where(less(warray,Pthresh)+greater(abs(Barray-Bideal),Dthresh),\
                            Bideal,Barray)
        #use a weighted least squares to place images (large weights for good correlation results, low weights for bad)
        posfit,resids = weighted_linear_least_squares(Aarray,Barray,warray**2) 
        #finally, store positions
        for j in range(len(othinds)):
            positions[cz_inds[othinds[j]],1]=posfit[2*j]
            positions[cz_inds[othinds[j]],2]=posfit[2*j+1]
    return positions

def run_CCimages(img_list,search_width=200,search_height=200,xoff=0.0,yoff=0.0,templ_width=50,templ_height=50):
    CCimages_results = gen_tempfile("CCimages_results","txt")
    cmdstr = "/home/bjnieman/source/TV/CCimages -w %d -l %d -x %f -y %f -t %d -u %d %s %s"%\
              (search_width,search_height,xoff,yoff,templ_width,templ_height,img_list[0],img_list[1])
    cmdstr = cmdstr + " > %s"%CCimages_results
    print cmdstr
    os.system(cmdstr)
    fh=open(CCimages_results)
    xoffnew=xoff; yoffnew=yoff;
    CCresult=0.0; Imean=0.0;
    while (1):
        currtext = fh.readline()
        if currtext=='': break
        print currtext
        reobj=re.search('offset_x=[0-9-]+',currtext)
        if not (reobj==None):
            xoffnew=int(reobj.group(0)[9:])
        reobj=re.search('offset_y=[0-9-]+',currtext)
        if not (reobj==None):
            yoffnew=int(reobj.group(0)[9:])
        reobj=re.search('CC=[0-9.]+',currtext)
        if not (reobj==None):
            CCresult=float(reobj.group(0)[4:])
        reobj=re.search('Imean=[0-9.]+',currtext)
        if not (reobj==None):
            Imean=float(reobj.group(0)[6:])
    fh.close()
    rmfilelist([CCimages_results])
    return array((xoffnew,yoffnew),float),CCresult,Imean

def run_image_overlay(imglist,positions,outimg_size_x=None,outimg_size_y=None):
    if (outimg_size_x==None):
        (matX,matY)=Image.open(imglist[0]).size
        min_offset_x = minimum.reduce(positions[:,-1])
        positions[:,-1] -= min_offset_x
        max_offset_x = maximum.reduce(positions[:,-1])
        outimg_size_x = max_offset_x + matX
    if (outimg_size_y==None):
        (matX,matY)=Image.open(imglist[0]).size
        min_offset_y = minimum.reduce(positions[:,-2])
        offset_results[:,-2] -= min_offset_y
        max_offset_y = maximum.reduce(positions[:,-2])
        outimg_size_y = max_offset_y + matY
    image_overlay_output = gen_tempfile("image_overlay_Z%04d"%positions[0,0],imglist[0][-3:])
    cmdstr = "/home/bjnieman/source/TV/image_overlay -x %d -y %d "%(int(outimg_size_x),int(outimg_size_y))
    for j in range(len(imglist)):
        cmdstr+="%s %f %f "%(imglist[j],positions[j,-1],positions[j,-2])
    cmdstr+=image_overlay_output
    print cmdstr
    os.system(cmdstr)
    return image_overlay_output

def save_positions_to_file(img_list,coord_array,positions,outputfile):
    print "Outputting %s...\n"%outputfile
    fh=open(outputfile,'w')
    for j in range(len(img_list)):
        fh.write("%s (%d %d %d) (%f %f %f)\n"%(img_list[j],\
                  coord_array[j,0],coord_array[j,1],coord_array[j,2],\
                  positions[j,0],positions[j,1],positions[j,2]))
    fh.close()
    return 1

def get_positions_from_file(positions_file):
    print "Reading positions from file (%s)...\n"%positions_file
    fh=open(positions_file,'r')
    entries=fh.readlines()
    fh.close()
    positions=[]
    coordlist=[]
    for cline in entries:
        reoutput = re.findall(r'[a-zA-Z0-9/\-_\.]+', cline)
        if (len(reoutput)==0): continue
        positions.append([reoutput[-3],reoutput[-2],reoutput[-1]])
        coordlist.append([reoutput[-6],reoutput[-5],reoutput[-4]])
    return array(positions,float), array(coordlist,int)

def generate_mnc_file_from_tifstack(Zstacklist,outputfile,zstep=0.01,ystep=0.00137,xstep=0.00137):
    outputprefix=outputfile[:-4]
    (newmatX,newmatY)=Image.open(Zstacklist[0]).size
    Zmnclist=[]; Graylist=[]
    j=0
    for cfile in Zstacklist:
        cgrayfile = gen_tempfile(outputprefix.split('/')[-1]+'_Z%04d'%j,'gray'); 
        cmncfile = gen_tempfile(outputprefix.split('/')[-1]+'_Z%04d'%j,'mnc');
        j+=1
        Graylist.append(cgrayfile); Zmnclist.append(cmncfile)
        cmdstr="convert -type Grayscale -size %dx%d %s %s"%(newmatY,newmatX,cfile,cgrayfile)
        print cmdstr
        os.system(cmdstr)
        cmdstr="cat %s | rawtominc -byte -unsigned -xstep %f -ystep %f -zstep %f -origin 0 0 0 %s 1 %d %d"%\
                (cgrayfile,float(xstep),float(ystep),float(zstep),cmncfile,newmatY,newmatX)
        print cmdstr
        os.system(cmdstr)
    cmdstr="mincconcat -clobber -2 -sequential "
    for cfile in Zmnclist: cmdstr=cmdstr+" "+cfile
    cmdstr=cmdstr+" "+outputfile
    print cmdstr
    os.system(cmdstr)
    cmdstr="minc_modify_header -dinsert zspace:step=0.01 %s"%outputfile
    print cmdstr
    os.system(cmdstr)
    rmfilelist(Zmnclist)
    rmfilelist(Graylist)
    return outputfile

def parse_mosaic_logfile(logfilename):
    fh=open(logfilename)
    while (1):
        currtext = fh.readline()
        if currtext=='': break
        reobj=re.search('rows:[0-9-]+',currtext)
        if not (reobj==None): nrows=int(reobj.group(0)[5:])
        reobj=re.search('columns:[0-9-]+',currtext)
        if not (reobj==None): ncolumns=int(reobj.group(0)[8:])
        reobj=re.search('sectionres:[0-9-]+',currtext)
        if not (reobj==None): sectionres=int(reobj.group(0)[11:])
    fh.close()
    xres = 0.001*[TV_LORES,TV_HIRES][nrows>832]
    yres = 0.001*[TV_LORES,TV_HIRES][ncolumns>832]
    zres = 0.001*sectionres  #this really needs to incorporate both zres and sectionres
    return xres,yres,zres

def rmfilelist(filelist):
    for junkfile in filelist:
        cmdstr = 'rm %s'%junkfile
        os.system(cmdstr)
    return None


if __name__ == '__main__':

    usage = """%s <input_directory> <output_mnc_file>
   or  %s --help

%s is a script that stitches together tiles from TV. Images are assumed
to lie on a grid as defined by the filenames (*Z#_Y#_X#*).
"""
    usage = usage % ((program_name, )*3)

    parser = OptionParser(usage)
    parser.add_option("--clobber", action="store_true", dest="clobber",
                       default=0, help="overwrite output file")
    parser.add_option("--use_positions_file", type="string", dest="use_positions_file", metavar="positions_file.txt", \
                      help="use an existing positions file instead of generating a new one")
    parser.add_option("--save_positions_file", type="string", dest="save_positions_file", metavar="positions_file.txt", \
                      help="save the final positions to file (for subsequent use with use_positions_file)")
    parser.add_option("--overlapx",type="float",dest="overlapx",default=20.0,\
                      help="tile overlap in percent")
    parser.add_option("--overlapy",type="float",dest="overlapy",default=15.0,\
                      help="tile overlap in percent")
    parser.add_option("--offsetguessy",type="float",dest="offsetguessy",default=-17.0,\
                      help="grid offset guess in y")
    parser.add_option("--offsetguessx",type="float",dest="offsetguessx",default=10.0,\
                      help="grid offset guess in x")
    parser.add_option("--input_prefix", type="string", dest="input_prefix", default="Tile", \
                      help="input prefix (default: 'Tile')")
    parser.add_option("--Zstart",type="int",dest="Zstart",default=-1,\
                      help="Z start index")
    parser.add_option("--Ystart",type="int",dest="Ystart",default=-1,\
                      help="Y start index")
    parser.add_option("--Xstart",type="int",dest="Xstart",default=-1,\
                      help="X start index")
    parser.add_option("--Zend",type="int",dest="Zend",default=-1,\
                      help="Z end index")
    parser.add_option("--Yend",type="int",dest="Yend",default=-1,\
                      help="Y end index")
    parser.add_option("--Xend",type="int",dest="Xend",default=-1,\
                      help="X end index")
    parser.add_option("--Zref",type="int",dest="Zref",default=-1,\
                      help="Z plane reference")
    parser.add_option("--Zsep",type="float",dest="Zsep",default=0.05,\
                      help="Z slice separation in mm")
    parser.add_option("--index_length",type="int",dest="index_length",default=3,\
                      help="number of characters in file index (default: 3, i.e. Z001)")
    parser.add_option("--file_type", type="string", dest="file_type", metavar="file_extension", \
                      default="tif",help="file extension (default: tif)")
    parser.add_option("--keepstack", action="store_true", dest="keepstack",
                       default=0, help="keep stack of images instead of generating mnc file (output name used as prefix)")
    parser.add_option("--Zstack_pzIcorr", action="store_true", dest="Zstack_pzIcorr",
                       default=0, help="intensity correct piezo stacked images")
    parser.add_option("--mosaic_logfile", type="string", dest="mosaic_logfile", metavar="mosaic_logfile.txt", \
                       default=None, help="TissueVission log file with acquisition info (for extracting resolution, etc.)")


    options, args = parser.parse_args()

    try:
        if len(args)<2: raise FatalError, "Specify input directory and output file."
        inputdirectory = args[-2]
        outputfile = args[-1]
        outputprefix = outputfile[:-4]
        if not options.clobber and os.path.exists(outputfile):
            raise FatalError, \
               "The --clobber option is needed to overwrite an existing file."
    except FatalError, e:
        print 'Error(%s):' % program_name, e.msg
        raise SystemExit

    #generate image list
    starts=[]
    for j in [options.Zstart,options.Ystart,options.Xstart]:
        if (j>=0): starts.append(j)
        else: starts.append(None)
    ends=[]
    for j in [options.Zend,options.Yend,options.Xend]:
        if (j>=0): ends.append(j)
        else: ends.append(None)
    parsed_img_list,coord_array=compose_img_file_list(inputdirectory,options.input_prefix,
                                                      starts=starts,
                                                      ends=ends,
                                                      index_length=options.index_length,
                                                      filetype=options.file_type)
    uniqueZ=unique(coord_array[:,0])
    (matX,matY)=Image.open(parsed_img_list[0]).size

    #determine offsets with CCimages or read in positions from previously written file
    existing_positions_file_flag = getattr(options,'use_positions_file')
    if not existing_positions_file_flag:
        positions=compute_offsets(parsed_img_list,coord_array,overlapx=options.overlapx,overlapy=options.overlapy,\
                                            offsety=options.offsetguessy,offsetx=options.offsetguessx,Zref=options.Zref)
        if getattr(options,'save_positions_file'):
            save_positions_to_file(parsed_img_list,coord_array,positions,options.save_positions_file)
    else:
        positions,file_coord_array=get_positions_from_file(options.use_positions_file)
        #we probably need to match file_coord_array with coord_array and sort the positions here

    #adjust positions to be all positive offsets based on global minima
    min_offset_x = minimum.reduce(positions[:,-1])
    min_offset_y = minimum.reduce(positions[:,-2])
    positions[:,-1] -= min_offset_x
    positions[:,-2] -= min_offset_y
    max_offset_x = maximum.reduce(positions[:,-1])
    outimg_size_x = max_offset_x + matX
    max_offset_y = maximum.reduce(positions[:,-2])
    outimg_size_y = max_offset_y + matY

    #overlay images with opencv based run_image_overlay (enblend was way too slow)        
    Zstacklist=[]
    for z in uniqueZ:
        zinds=nonzero(coord_array[:,0]==z)[0]
        clist=[parsed_img_list[j] for j in zinds]
        #generate full slices from tiles
        Zsliceimg=run_image_overlay(clist,positions[zinds,:],outimg_size_x=outimg_size_x,outimg_size_y=outimg_size_y)
        Zstacklist.append(Zsliceimg)

    #if needed, perform slice-by-slice intensity normalization (i.e., for piezo stacks)
    #NEED TO SORT AND PROCESS STACK IN PIEZO BATCHES
    if (options.Zstack_pzIcorr):
        intensity_normalize_Zstack(Zstacklist,Zstacklist)

    #generate minc file
    #below: need to read in step sizes from somewhere? a TissueVision log file?
    if (options.keepstack):
        for z,cfile in enumerate(Zstacklist):
            cmdstr="cp %s %s"%(cfile,outputfile+'_Z%04d'%uniqueZ[z]+'.%s'%options.file_type)
            os.system(cmdstr)
    else:
        (x_step,y_step,z_step) = (0.00137,0.00137,options.Zsep) #default resolutions
        if (options.mosaic_logfile!=None):
            x_step,y_step,z_step = parse_mosaic_logfile(options.mosaic_logfile)
        generate_mnc_file_from_tifstack(Zstacklist,outputfile,zstep=z_step,ystep=y_step,xstep=x_step) 

    #clean up temp files
    rmfilelist(Zstacklist)

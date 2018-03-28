#!/usr/bin/env python3
#
# TV_stitch.py
#
# stitch XY planes of TV tiles
#
# Created February 2012

from sys import argv,path
import string
import os
import subprocess
import getopt
from optparse import OptionParser, Option, OptionValueError
import re
from numpy import *
from numpy.linalg import lstsq
import glob
import operator

from python.Zstack_icorr import *

from scipy.stats import t

program_name = 'TV_stitch.py'

#CURRENTLY TOGGLES BETWEEN ONLY A LO AND A HI RES MODE
#NEED TO MAKE THIS MORE FLEXIBLE
TV_LORES = 1.37  #um, calibrated
TV_HIRES = 0.55  #um, calibrated
SHAVE_LORES = 15 #pixels to crop off edges
SHAVE_HIRES = 38 
LORESMAT = 832  
HIRESMAT = 2080
STAGE_CALIB_FACTOR = 0.1 #um (per step)

TEMPDIRECTORY = "./tmp_"+program_name+"_"+str(os.getpid()) #"/tmp/"+program_name+'_'+str(os.getpid())
VERBOSE = False

MIN_WEIGHT = 1e-3

#----------------------------------------------------------------------------
# define program specific exception
class FatalError(BaseException):
    def __init__(self,args=None):
        self.msg = args

# define class with relevant tile info (indexed position, reported position)
class Tile(object):
    def __init__(self, filename=None, croppedfilename=None, croppedfilteredfilename=None, processedfilename=None, \
                       indexarray=[-1,-1,-1], posarray=[nan,nan,nan], pixoffsetarray=[nan,nan,nan]):
        self.filename = filename; self.croppedfilename = croppedfilename; self.croppedfilteredfilename=croppedfilteredfilename; self.processedfilename=processedfilename
        self.indexarray = array(indexarray,int); self.posarray = array(posarray,float)
        self.pixoffsetarray = array(pixoffsetarray,float)
    def __str__(self):
        return "filename = %s\ncroppedfilename = %s\nindexarray = [%d, %d, %d]\nposarray = [%f, %f, %f]\npixoffsetarray = [%f, %f, %f]\n" % \
               (self.filename,self.croppedfilename,self.indexarray[0],self.indexarray[1],self.indexarray[2], \
                self.posarray[0],self.posarray[1],self.posarray[2],self.pixoffsetarray[0],self.pixoffsetarray[1],self.pixoffsetarray[2])
    def __repr__(self):
        return "filename = %s\ncroppedfilename = %s\nindexarray = [%d, %d, %d]\nposarray = [%f, %f, %f]\npixoffsetarray = [%f, %f, %f]\n" % \
               (self.filename,self.croppedfilename,self.indexarray[0],self.indexarray[1],self.indexarray[2], \
                self.posarray[0],self.posarray[1],self.posarray[2],self.pixoffsetarray[0],self.pixoffsetarray[1],self.pixoffsetarray[2])

#---------------------------------------------------------------------------

def gen_tempfile(descrip_str,ftype_str):
    if not (os.path.exists(TEMPDIRECTORY)): os.mkdir(TEMPDIRECTORY)
    tempstr = TEMPDIRECTORY + "/" + program_name + '_' + descrip_str + '.' + ftype_str
    return tempstr

def run_subprocess(cmdstr):
    if VERBOSE:
        print(cmdstr)
    p=subprocess.Popen(cmdstr,stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    p.wait()
    (out,err)=p.communicate()
    if not (err==''):
        print(err)
    if VERBOSE:
        print(out)
    return out

def crop_img(infile,outfile,shavewidth=None,shaveheight=None,imgres='LORES',depth=None,im=False):
    if (shavewidth==None): shavewidth=[SHAVE_LORES,SHAVE_HIRES][ {'LORES':0, 'HIRES':1}[imgres] ]
    if (shaveheight==None): shaveheight=[SHAVE_LORES,SHAVE_HIRES][ {'LORES':0, 'HIRES':1}[imgres] ]
    if (im):
        depthstr = "-depth %d"%depth if (depth!=None) else ""
        cmdstr="convert -shave %dx%d %s %s %s"%(shavewidth,shaveheight,depthstr,infile,outfile)
    else:
        cmdstr="cvRectCrop -l %d -r %d -t %d -b %d %s %s"%(shavewidth,shavewidth,shaveheight,shaveheight,infile,outfile)
    cmdout = run_subprocess(cmdstr)
    return 0

def gen_avg_of_tiles(globstr,outfile,maxfiles=20000,Iscale=50.0):
    cmdstr='cvAvgImages -g "%s" -s %f -m %d %s'%(globstr,Iscale,maxfiles,outfile)
    cmdout = run_subprocess(cmdstr)
    return 0

def corr_tiles(globstr,avgTileImg,median_kernel_size=5,postfix="_avgIcorr"):
    cmdstr='cvCorrTiles -g "%s" -a %s -p %s -m %d'%(globstr,avgTileImg,postfix,median_kernel_size)
    cmdout = run_subprocess(cmdstr)
    return 0

def TV_parameters(paramfile):
    fh = open(paramfile,'r')
    text_lines = fh.readlines()
    fh.close()
    param_dict={}
    XPos=[]
    YPos=[]
    for curr_line in text_lines:
        m = re.search('^[XY]Pos',curr_line)
        if (m):
            if (curr_line[0]=='X'):
                XPos.append(int(curr_line.partition(':')[-1]))
            elif (curr_line[0]=='Y'):
                YPos.append(int(curr_line.partition(':')[-1]))
        else:
            words = curr_line.partition(':')
            varname = words[0]
            varvalue = words[-1].lstrip().rstrip()
            try:
                param_dict[varname]=int(varvalue)
            except ValueError:
                try:
                    param_dict[varname]=float(varvalue)
                except ValueError:
                    param_dict[varname]=varvalue
    posarray = empty((len(XPos),2),int)
    posarray[:,0] = array(XPos)
    posarray[:,1] = array(YPos)
    param_dict['posarray'] = STAGE_CALIB_FACTOR * posarray #gives reported stage position in um
    return param_dict

numbers = re.compile(r'([0-9.]+)')
def numericalSort(value):
    parts = numbers.split(value)
    fileintvals = list(map(float, parts[1::2]))
    return fileintvals[-2]

#def generate_gradcombined_images(infile,outfile,combineflag=True,im=False):
#    cfile=infile
#    if (combineflag):
#        cgradfile = gen_tempfile(infile.split('/')[-1][:-4]+'_grad',outfile[-3:])
#    else:
#        cgradfile = outfile
#    if (im):
#        cmdstr = "convert %s -define convolve:scale=50%%^ -bias 50%% -convolve '9x9: 0,1,1,2,2,2,1,1,0, 1,2,4,5,5,5,4,2,1, 1,4,5,3,0,3,5,4,1, 2,5,3,-12,-24,-12,3,5,2, 2,5,0,-24,-40,-24,0,5,2, 2,5,3,-12,-24,-12,3,5,2, 1,4,5,3,0,3,5,4,1, 1,2,4,5,5,5,4,2,1, 0,1,1,2,2,2,1,1,0' %s"%(cfile,cgradfile)
#    else:
#        cmdstr = "cvFilter -k scharr -w 9 -s 1.4 %s %s"%(cfile,cgradfile)
#    cmdout = run_subprocess(cmdstr)
#    if (combineflag):
#        if (im):
#            cmdstr = "convert %s %s -set colorspace RGB -combine %s"%(cfile,cgradfile,outfile)
#        else:
#            cmdstr = "cvMerge -r %s -g %s %s"%(cfile,cgradfile,outfile)
#        cmdout = run_subprocess(cmdstr)
#    return 0

def generate_gradcombined_images(infile,outfile,combineflag=True,im=False):
    cfile=infile
    if (im): 
        if (combineflag):
            cgradfile = gen_tempfile(infile.split('/')[-1][:-4]+'_grad',outfile[-3:])
        else:
            cgradfile = outfile
            cmdstr = "convert %s -define convolve:scale=50%%^ -bias 50%% -convolve '9x9: 0,1,1,2,2,2,1,1,0, 1,2,4,5,5,5,4,2,1, 1,4,5,3,0,3,5,4,1, 2,5,3,-12,-24,-12,3,5,2, 2,5,0,-24,-40,-24,0,5,2, 2,5,3,-12,-24,-12,3,5,2, 1,4,5,3,0,3,5,4,1, 1,2,4,5,5,5,4,2,1, 0,1,1,2,2,2,1,1,0' %s"%(cfile,cgradfile)
        cmdout = run_subprocess(cmdstr)
        if (combineflag):
            cmdstr = "convert %s %s -set colorspace RGB -combine %s"%(cfile,cgradfile,outfile)
            cmdout = run_subprocess(cmdstr)
    else:
        cgradfilex = gen_tempfile(infile.split('/')[-1][:-4]+'_gradx',outfile[-3:])
        cgradfiley = gen_tempfile(infile.split('/')[-1][:-4]+'_grady',outfile[-3:])
        cmdstr = "cvFilter -k scharrx -w 9 -s 1.4 %s %s"%(cfile,cgradfilex)
        cmdout = run_subprocess(cmdstr)
        cmdstr = "cvFilter -k scharry -w 9 -s 1.4 %s %s"%(cfile,cgradfiley)
        cmdout = run_subprocess(cmdstr)
        if (combineflag):
            cmdstr = "cvMerge -r %s -g %s -b %s %s"%(cfile,cgradfilex,cgradfiley,outfile)
        else:
            cmdstr = "cvMerge -g %s -b %s %s"%(cgradfilex,cgradfiley,outfile)
        cmdout = run_subprocess(cmdstr)
    return 0

def run_cvFilter(infile,outfile,kernelstr="gauss",kernelwidth=9,filtwidth=1.4):
    cmdstr = "cvFilter -k %s -w %i -s %f %s %s"%(kernelstr,kernelwidth,filtwidth,infile,outfile)
    cmdout = run_subprocess(cmdstr)
    return 0

def generate_preprocessed_images(inputdirectory,starts=[None,None,None],ends=[None,None,None],channelflag=1,\
                                 imgftype='tif',fastpiezoloop=False,gradcombine=False,im=False,
                                 corr_tile_nonuniformity=False,medfilter_tile=False,medfilter_size=3):
    #find and compose list of directories to work with
    try:
        fulldirectorylist = glob.glob(inputdirectory+'-[0-9]*')
        if not fulldirectorylist:
            raise FatalError("Cannot find input directory(ies).")
    except FatalError as e:
        print('Error(%s):' % program_name, e.msg)
        raise SystemExit
    #refine and sort directory list
    direclist=[]
    inputdirhead,junk,input_prefix = inputdirectory.rpartition('/')
    reprog=re.compile(input_prefix+'[0-9-]*\Z')
    for cname in fulldirectorylist:
        m=reprog.search(cname)
        if (m!=None): direclist.append(inputdirhead+'/'+m.group(0))
    direclist.sort()
    #read #X,#Y,#Zpeizo from first Mosaic text file
    TVscanfile = direclist[0]+'/'+'Mosaic_'+direclist[0].rsplit('/')[-1]+'.txt'
    TVparamdict = TV_parameters(TVscanfile)
    (N_x,N_y,N_z_slices,N_z_piezo) = (TVparamdict['mcolumns'],TVparamdict['mrows'],\
                                      TVparamdict['sections'],[1,TVparamdict['layers']][TVparamdict['Zscan']])
    if (len(direclist)!=N_z_slices):
        print("Mismatch in Mosaic file 'sections' (%d) and number of identified directories (%d)"%(N_z_slices,len(direclist)))
        N_z_slices=len(direclist)
        TVparamdict['sections']=N_z_slices
    try:
        globlist = glob.glob(direclist[0]+'/'+'*-*-*'+'_%02d.'%channelflag+imgftype)
        if (len(globlist)==0):
            raise FatalError("Cannot find files for the specified channel.")
    except FatalError as e:
        print('Error(%s):' % program_name, e.msg)
        raise SystemExit
    #due to weird behaviour of convert in crop_img for some images, force depth to match input (so it doesn't change)
    imgdepth = int( run_subprocess("identify -format \"%%z\" %s"%globlist[0]) )
    #generate complete file list with indexed and reported positions
    TileList=[]
    for j,cname in enumerate(direclist):
        #read Mosaic file for current directory to get position info
        TVscanfile = cname+'/'+'Mosaic_'+cname.rsplit('/')[-1]+'.txt'
        TVparamdict = TV_parameters(TVscanfile)
        globlist = glob.glob(cname+'/'+'*-*-*'+'_%02d.'%channelflag+imgftype)
        globlist.sort(key=numericalSort)
        nfiles = len(globlist)
        if (nfiles!=N_x*N_y*N_z_piezo):
            print("******************************************************************")
            print("WARNING: Incorrect number of files identified in %s"%cname)
            print("******************************************************************")
            if (nfiles>N_x*N_y*N_z_piezo): nfiles=N_x*N_y*N_z_piezo
        if (nfiles/N_z_piezo!=TVparamdict['posarray'].shape[0]):
            print("**************************************************************************")
            print("WARNING: Inconsistent number of files (%d) for Mosaic position list (%d)"%(len(globlist),TVparamdict['posarray'].shape[0]))
            print("**************************************************************************")
            nfiles = min([nfiles,N_z_piezo*TVparamdict['posarray'].shape[0]])
        globlist=globlist[0:nfiles] #this line here only because sometimes TV errors generate extra files
        for k in range(nfiles):
            cfile = globlist[k]
            if (fastpiezoloop):
                zindex = j*N_z_piezo+1+[0,k%N_z_piezo][N_z_piezo>1]
                posindex = k/N_z_piezo
                try:
                    reported_zpos = 0.001*TVparamdict['sectionres']*j + 0.001*2.0*TVparamdict['zres']*[0,k%N_z_piezo][N_z_piezo>1] #2X?
                except IndexError:
                    print(k, cfile)
                    raise SystemExit
            else:
                zindex = j*N_z_piezo+1+[0,k/(N_x*N_y)][N_z_piezo>1]
                posindex = k%(N_x*N_y)
                reported_zpos = 0.001*TVparamdict['sectionres']*j + 0.001*2.0*TVparamdict['zres']*[0,k/(N_x*N_y)][N_z_piezo>1] #2X?       
            reported_xpos = TVparamdict['posarray'][posindex,1]    #reported position from Mosaic file
            reported_ypos = TVparamdict['posarray'][posindex,0]*-1.0    #but note swap of y-->-x and x-->y due to TV versus image convention
            TileList.append( Tile(cfile,None,None,None,[zindex,0,0],[reported_zpos,reported_ypos,reported_xpos]) )
    #assign indices, sort files
    xmin = min( [ctile.posarray[2] for ctile in TileList] ); xmax = max( [ctile.posarray[2] for ctile in TileList] )
    ymin = min( [ctile.posarray[1] for ctile in TileList] ); ymax = max( [ctile.posarray[1] for ctile in TileList] )
    for ctile in TileList:
        ctile.indexarray[2] = 1+round( (N_x-1)*(ctile.posarray[2]-xmin)/float(xmax-xmin) )
        ctile.indexarray[1] = 1+round( (N_y-1)*(ctile.posarray[1]-ymin)/float(ymax-ymin) )
    f_indextuple = lambda ctile: (ctile.indexarray[0],ctile.indexarray[1],ctile.indexarray[2])
    TileList.sort( key=f_indextuple )
    #now crop images, working only within specified start and end
    for j in range(len(TileList)):
        #skip to next image if outside desired start:end range
        for k in range(3):
            if (starts[k]!=None): #'<' not supported between 'int' and 'NoneType'
                if (TileList[j].indexarray[k]<starts[k]): continue
            if (ends[k]!=None):
                if (TileList[j].indexarray[k]>ends[k]): continue

        #crop
        TileList[j].croppedfilename = gen_tempfile('Tile_Z%03d_Y%03d_X%03d'%(TileList[j].indexarray[0],\
                              TileList[j].indexarray[1],TileList[j].indexarray[2]),TileList[j].filename.split('.')[-1])
        if not (os.path.exists(TileList[j].croppedfilename)):
            crop_img(TileList[j].filename,TileList[j].croppedfilename,imgres=['LORES','HIRES'][TVparamdict['rows']>LORESMAT],depth=imgdepth,im=im)
    if (corr_tile_nonuniformity):
       avgTileImg = gen_tempfile('avgTileImg',TileList[0].filename.split('.')[-1])
       globstr = os.path.join(TEMPDIRECTORY,program_name+"_Tile_Z[0-9][0-9][0-9]_Y[0-9][0-9][0-9]_X[0-9][0-9][0-9]."+TileList[0].filename.split('.')[-1])
       if not (os.path.exists(avgTileImg)):
           gen_avg_of_tiles(globstr,avgTileImg,maxfiles=20000,Iscale=50.0)
           postfix="_avgIcorr"
           corr_tiles(globstr,avgTileImg,median_kernel_size=5,postfix=postfix)
           for j in range(len(TileList)):
               if not (TileList[j].croppedfilename is None):
                   TileList[j].croppedfilename = TileList[j].croppedfilename[:-4]+postfix+'.'+TileList[j].filename.split('.')[-1]
    for j in range(len(TileList)):
        if (medfilter_tile):
            if not (TileList[j].croppedfilename is None):
                TileList[j].croppedfilteredfilename = TileList[j].croppedfilename[:-4]+"_medfilt"+'.'+TileList[j].filename.split('.')[-1]
                #run median filter
                run_cvFilter(TileList[j].croppedfilename,TileList[j].croppedfilteredfilename,kernelstr="median",kernelwidth=medfilter_size)
        else:
            TileList[j].croppedfilteredfilename = TileList[j].croppedfilename
    for j in range(len(TileList)):
        if (TileList[j].croppedfilename is None): 
            continue
        if (gradcombine):
            TileList[j].processedfilename = gen_tempfile('Tile_gradRGB_Z%03d_Y%03d_X%03d'%(TileList[j].indexarray[0],\
                              TileList[j].indexarray[1],TileList[j].indexarray[2]),TileList[j].filename.split('.')[-1])
            if not (os.path.exists(TileList[j].processedfilename)):
                generate_gradcombined_images(TileList[j].croppedfilteredfilename,TileList[j].processedfilename,combineflag=True,im=im)
        else:
            TileList[j].processedfilename = TileList[j].croppedfilteredfilename
        #pixoffset guess based on only TV Mosaic positions
        TileList[j].pixoffsetarray[2] = TileList[j].posarray[2]/[TV_LORES,TV_HIRES][TVparamdict['rows']>LORESMAT]
        TileList[j].pixoffsetarray[1] = TileList[j].posarray[1]/[TV_LORES,TV_HIRES][TVparamdict['rows']>LORESMAT]
        TileList[j].pixoffsetarray[0] = TileList[j].indexarray[0]
    #delete TileList elements with images that were not cropped
    TileList = [ctile for ctile in TileList if ctile.croppedfilename!=None]
    TVparamdict['N_z_piezo'] = N_z_piezo
    return TileList,TVparamdict

def weighted_linear_least_squares(Amat,bcol,w):
    matrixmultiply=dot
    npts = len(Amat)
    W = identity(npts,dtype=float)
    for j in range(npts):
        W[j,j]=w[j]
    Bmat= bcol
    WAmat = matrixmultiply(W,Amat)
    Wbcol = matrixmultiply(W,bcol)
    x,yresids,rank,s = lstsq(WAmat,Wbcol)
    return x,yresids

def compute_offsets(TileList,overlapx=20.0,overlapy=15.0,Zref=-1,Cthresh=0.3,zsearch=40):
    #determine offsets of a grid of images with overlap at the edges based on CCimages (opencv correlative template matching)
    Nfiles=len(TileList)
    cmdout=run_subprocess("identify -format \"%%w %%h\" %s"%TileList[0].croppedfilename)
    matX=int(cmdout.split()[0]); matY=int(cmdout.split()[1])
    uniqueZ = unique([ctile.indexarray[0] for ctile in TileList])
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
    templx_dict={'x': 0.5*matX*overlapx/100.0,  'y': matY-2*searchy_dict['x'], 'z': matX-2*zsearch}
    temply_dict={'x': matX-2*searchx_dict['y'], 'y': 0.5*matY*overlapy/100.0,  'z': matY-2*zsearch}
    areamax = float( max([templx_dict['x']*temply_dict['x'],templx_dict['y']*temply_dict['y'],templx_dict['z']*temply_dict['z']]) )
    axis_weights={'x': sqrt(templx_dict['x']*temply_dict['x']/areamax), 
                  'y': sqrt(templx_dict['y']*temply_dict['y']/areamax), 
                  'z': sqrt(templx_dict['z']*temply_dict['z']/areamax) } #weight correlations based on area of template
    grid_offset_results=array([[0,0,0,0]],float)
    refPmatch_results=[]
    for zindex,currz in enumerate(sorted_uniqueZ):
        cz_inds = [cind for cind,ctile in enumerate(TileList) if ctile.indexarray[0]==currz]
        cz_nfiles = len(cz_inds)
        Alist=[]     #coefficients of positions to produce offsets (i.e. 0 0 0 ... 1 -1 0 ... 0)
        blist=[]     #offset results
        wlist=[]     #stores CCimages correlation results, serves as weights for least squares fit
        Ilist=[]     #stores Imean for CCimages results
        zposlist=[]  #the x or y position from the previous z-slice (NOT the z-position)
        direclist=[] #string with direction and coordinates
        for caxis in ['x','y','z']:
            searchx=int(searchx_dict[caxis])
            searchy=int(searchy_dict[caxis])
            templx =int(templx_dict[caxis])
            temply =int(temply_dict[caxis])
            for j in range(0,cz_nfiles):
                cfile=TileList[cz_inds[j]].croppedfilename
                #identify current reference image
                relpos = array( [(TileList[cz_inds[j]].indexarray - ctile.indexarray) for ctile in TileList] )
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
                #run CCimages
                if (caxis=='z'):
                    (xoff,yoff) = (0.0,0.0)
                else:
                    xoff = TileList[cz_inds[j]].pixoffsetarray[2]-TileList[refind].pixoffsetarray[2]
                    yoff = TileList[cz_inds[j]].pixoffsetarray[1]-TileList[refind].pixoffsetarray[1]
                newoffsets,CCresult,Imean=run_CCimages([TileList[refind].processedfilename,TileList[cz_inds[j]].processedfilename],\
                                                                       search_width=searchx,\
                                                                       search_height=searchy,\
                                                                       xoff=xoff,\
                                                                       yoff=yoff,\
                                                                       templ_width=templx,\
                                                                       templ_height=temply)
                #put results into Alist, blist, wlist and zposlist
                Arow=zeros((2*cz_nfiles,),float)
                Arow[2*j]=1.0 
                if not (caxis=='z'):
                    Arow[2*where(cz_inds==refind)[0][0]]=-1.0
                else:
                    newoffsets[1]+=TileList[refind].pixoffsetarray[1]
                zposlist.append(TileList[refind].pixoffsetarray[1])
                direclist.append(caxis+'y')
                Alist.append(Arow); blist.append(newoffsets[1])
                wlist.append(CCresult*axis_weights[caxis]); Ilist.append(Imean) 
                Arow=zeros((2*cz_nfiles,),float)
                Arow[2*j+1]=1.0
                if not (caxis=='z'):
                    Arow[2*where(cz_inds==refind)[0][0]+1]=-1.0
                else:
                    newoffsets[0]+=TileList[refind].pixoffsetarray[2]
                zposlist.append(TileList[refind].pixoffsetarray[2])
                direclist.append(caxis+'x')
                Alist.append(Arow); blist.append(newoffsets[0])
                wlist.append(CCresult*axis_weights[caxis]); Ilist.append(Imean)
        #QC on wlist and Ilist
        maxI = max(Ilist)
        #import matplotlib as mpl
        #mpl.use('Agg')
        #import matplotlib.pyplot as pl
        #pl.hist(Iarray,256)
        #pl.savefig("Iarrayhist.pdf")
        #wlist = [ [x,MIN_WEIGHT][x>0.98] for x in wlist] #zero tiles or spikes can produce CC=1.0, but aren't useful
        #generate arrays for lsq
        if (zindex==0): #for first slice, need to define a reference arbitrarily
            refind = 0  #TileList is already sorted based on indexarray, so Y001_X001 should be first in list
            Arow=zeros((2*cz_nfiles,),float); Arow[0]=1.0
            Alist.append(Arow); blist.append(TileList[cz_inds[0]].pixoffsetarray[1]) 
            wlist.append(1.0); Ilist.append(maxI)
            direclist.append('zy'); zposlist.append(TileList[cz_inds[0]].pixoffsetarray[1])
            Arow=zeros((2*cz_nfiles,),float); Arow[1]=1.0
            Alist.append(Arow); blist.append(TileList[cz_inds[0]].pixoffsetarray[2])
            wlist.append(1.0); Ilist.append(maxI)
            direclist.append('zx'); zposlist.append(TileList[cz_inds[0]].pixoffsetarray[2])
        Aarray = array(Alist,float) #ADD
        Barray = array(blist,float); 
        warray = array(wlist,float);
        warray = where(warray>MIN_WEIGHT, warray,MIN_WEIGHT)
        Iarray = array(Ilist)
        minI = min(Iarray)
        Ilimit = 0.015*(maxI-minI)+minI
        warray = where(Iarray<Ilimit,MIN_WEIGHT,warray)
        #warray_orig=warray.copy()
        if (cz_nfiles>=8):
            refPmatch_results.append(median(warray))            
            Pthresh=Cthresh*median(refPmatch_results)   
            medij={}; stdij={}; dfij={}
            for cdi in ['xx','yy','xy','yx']:
                ijlist=[bi for di,bi,wi in zip(direclist,blist,wlist) if ((di == cdi) and (wi > Pthresh))]
                medij[cdi] = mean(ijlist); stdij[cdi] = std(ijlist); dfij[cdi] = len(ijlist)-1
            for cdi in ['zx','zy']:
                ijlist=[(bi-zi) for di,bi,wi,zi in zip(direclist,blist,wlist,zposlist) if ((di == cdi) and (wi > Pthresh))]
                medij[cdi] = mean(ijlist); stdij[cdi] = std(ijlist); dfij[cdi] = len(ijlist)-1
            Bideal = array([{'xx':medij['xx'],'xy':medij['xy'],'yy':medij['yy'],'yx':medij['yx'],'zx':zi,'zy':zi}[di] \
                            for di,zi in zip(direclist,zposlist)],float)
            Bstd = array([stdij[di] for di in direclist],float)
            Bdf = array([dfij[di] for di in direclist],float)
            wdist = 2.0*(1.0-t.cdf(abs(Barray - Bideal)/Bstd,Bdf))
            wdist = where(isnan(wdist),1.0,wdist)
            warray = warray * wdist
            #Barray[:]=where(less(warray,Pthresh),Bideal,Barray)
            #warray[:]=where(less(warray,Pthresh),MIN_WEIGHT,warray)
            #inds = nonzero( greater(abs(Barray-Bideal),Dthresh) )[0]
            #inds = unique( append(inds,1-2*(inds%2)) ) #need to throw out both Y and X offsets for one bad distance measure
            #Barray[inds] = Bideal[inds]
            #warray[inds] = MIN_WEIGHT
        #use a weighted least squares to place images (large weights for good correlation results, low weights for bad)
        warray = where(less(warray,MIN_WEIGHT),MIN_WEIGHT,warray) 
        posfit,resids = weighted_linear_least_squares(Aarray,Barray,warray)
        if (isnan(posfit).any()):   #LSQ failed! This shouldn't happen
            print("LSQ fail (Z %d)..."%currz)
            Aarray.tofile(TEMPDIRECTORY+"/"+"LSQfail_Aarray_Z%d"%currz)
            Barray.tofile(TEMPDIRECTORY+"/"+"LSQfail_Barray_Z%d"%currz)
            warray.tofile(TEMPDIRECTORY+"/"+"LSQfail_warray_Z%d"%currz)
            raise SystemExit
        #finally, store positions
        for j in range(len(cz_inds)):
            TileList[cz_inds[j]].pixoffsetarray[1]=posfit[2*j]
            TileList[cz_inds[j]].pixoffsetarray[2]=posfit[2*j+1]
    return 0

def run_CCimages(img_list,search_width=200,search_height=200,xoff=0.0,yoff=0.0,templ_width=50,templ_height=50,CCvsback=0):
    cmdstr = "CCimages -w %d -l %d -x %f -y %f -t %d -u %d %s %s"%\
              (search_width,search_height,xoff,yoff,templ_width,templ_height,img_list[0],img_list[1])
    if (CCvsback):
        cmdstr += " -n"
    cmdout = run_subprocess(cmdstr)
    #assume output from CCimages of form: Imean=0.196748, offset_x=691, offset_y=5, CC=0.044111 
    numlist=re.findall(r'=[ 0-9-.]+',cmdout)  
    Imean=float(numlist[0][1:])
    xoffnew=int(numlist[1][1:])
    yoffnew=int(numlist[2][1:])
    CCresult=float(numlist[3][1:])
    return array((xoffnew,yoffnew),float),CCresult,Imean

def run_image_overlay(imglist,positions,outimg_size_x=None,outimg_size_y=None,outputfiletype=None,outscale=None):
    if (outimg_size_x==None):
        cmdout=run_subprocess("identify -format \"%%w %%h\" %s"%imglist[0])
        matX=int(cmdout.split()[0]); matY=int(cmdout.split()[1])
        min_offset_x = minimum.reduce(positions[:,-1])
        positions[:,-1] -= min_offset_x
        max_offset_x = maximum.reduce(positions[:,-1])
        outimg_size_x = max_offset_x + matX
    if (outimg_size_y==None):
        cmdout=run_subprocess("identify -format \"%%w %%h\" %s"%imglist[0])
        matX=int(cmdout.split()[0]); matY=int(cmdout.split()[1])
        min_offset_y = minimum.reduce(positions[:,-2])
        offset_results[:,-2] -= min_offset_y
        max_offset_y = maximum.reduce(positions[:,-2])
        outimg_size_y = max_offset_y + matY
    image_overlay_output = gen_tempfile("image_overlay_Z%04d"%positions[0,0],imglist[0][-3:])
    cmdstr = "image_overlay -x %d -y %d "%(int(outimg_size_x),int(outimg_size_y))
    if (outputfiletype=='byte'):
        cmdstr += "-s "
    elif (outputfiletype=='short'):
        cmdstr += "-i "
    if not (outscale==None):
        cmdstr += "-I %f "%outscale
    for j in range(len(imglist)):
        cmdstr+="%s %f %f "%(imglist[j],positions[j,-1],positions[j,-2])
    cmdstr+=image_overlay_output
    cmdout = run_subprocess(cmdstr)
    return image_overlay_output

def save_positions_to_file(TileList,outputfile):
    print("Outputting %s...\n"%outputfile)
    fh=open(outputfile,'w')
    for j in range(len(TileList)):
        fh.write("%s (%d %d %d) (%f %f %f)\n"%(TileList[j].filename,\
                  TileList[j].indexarray[0],TileList[j].indexarray[1],TileList[j].indexarray[2],\
                  TileList[j].pixoffsetarray[0],TileList[j].pixoffsetarray[1],TileList[j].pixoffsetarray[2])
                )
    fh.close()
    return 1

def get_positions_from_file(TileList,positions_file):
    print("Reading positions from file (%s)...\n"%positions_file)
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
    for ctile in TileList:
        matching_index = next((j for j,c_coord in enumerate(coordlist) if (array(c_coord,int)==ctile.indexarray).all()),-1)
        if (matching_index>=0): 
            ctile.pixoffsetarray = array(positions[matching_index],float)
        else:
            print("Failed to find matching coordinate indices for %s"%ctile.filename)
    return 0

def generate_mnc_file_from_tifstack(Zstacklist,outputfile,zstep=0.01,ystep=TV_LORES,xstep=TV_LORES,outdatatype="byte",Zcoordlist=None):
    outputprefix=outputfile[:-4]
    cmdout=run_subprocess("identify -format \"%%w %%h\" %s"%Zstacklist[0])
    newmatX=int(cmdout.split()[0]); newmatY=int(cmdout.split()[1])
    Zmnclist=[]; Graylist=[]
    if (Zcoordlist==None):
        Zcoordlist=zeros(len(Zstacklist),float)
        mncseqflag = "-sequential"
    else:
        mncseqflag = ""
    j=0
    for k,cfile in enumerate(Zstacklist):
        cgrayfile = gen_tempfile(outputprefix.split('/')[-1]+'_Z%04d'%j,'gray'); 
        cmncfile = gen_tempfile(outputprefix.split('/')[-1]+'_Z%04d'%j,'mnc');
        j+=1
        Graylist.append(cgrayfile); Zmnclist.append(cmncfile)
        cmdstr="convert -type Grayscale -size %dx%d %s %s"%(newmatY,newmatX,cfile,cgrayfile)
        cmdout = run_subprocess(cmdstr)
        cmdstr="cat %s | rawtominc -%s -unsigned -xstep %f -ystep %f -zstep %f -origin 0 0 %f %s 1 %d %d"%\
                (cgrayfile,outdatatype,float(xstep),float(ystep),float(zstep),Zcoordlist[k],cmncfile,newmatY,newmatX)
        cmdout = run_subprocess(cmdstr)
    cmdstr="mincconcat -clobber -2 %s "%mncseqflag
    for cfile in Zmnclist: cmdstr=cmdstr+" "+cfile
    cmdstr=cmdstr+" "+outputfile
    cmdout = run_subprocess(cmdstr)
    if (mncseqflag!=""):
        cmdstr="minc_modify_header -dinsert zspace:step=%f %s"%(zstep,outputfile)
        cmdout = run_subprocess(cmdstr)
    rmfilelist(Zmnclist)
    rmfilelist(Graylist)
    return outputfile

def rmfilelist(filelist):
    for junkfile in filelist:
        cmdstr = 'rm %s'%junkfile
        cmdout = run_subprocess(cmdstr)
    return None

####################################################################################################################################

if __name__ == '__main__':

    usage = """%s <input_directory> <output_mnc_file>
   or  %s --help

%s is a script that stitches together tiles from TV. Images are assumed
to come straight off the TissueVision system and undergo preprocessing.
"""
    usage = usage % ((program_name, )*3)

    parser = OptionParser(usage)
    parser.add_option("--clobber", action="store_true", dest="clobber",
                       default=False, help="overwrite output file")
    parser.add_option("--gradimag", action="store_true", dest="gradimag",
                       default=True, help="use gradient and raw image combined for correlation (default)")
    parser.add_option("--nogradimag", action="store_false", dest="gradimag",
                       default=True, help="do not use gradient and raw image combined for correlation")
    parser.add_option("--use_positions_file", type="string", dest="use_positions_file", metavar="positions_file.txt", \
                      help="use an existing positions file instead of generating positions from the input")
    parser.add_option("--save_positions_file", type="string", dest="save_positions_file", metavar="positions_file.txt", \
                      help="save the final positions to file (for subsequent use with use_positions_file)")
    parser.add_option("--overlapx",type="float",dest="overlapx",default=20.0,\
                      help="tile overlap in percent")
    parser.add_option("--overlapy",type="float",dest="overlapy",default=20.0,\
                      help="tile overlap in percent")
    parser.add_option("--channel",type="int",dest="channel",default=1,\
                      help="channel to stitch")
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
                      help="Z plane reference during tiling")
    parser.add_option("--Zstack_pzIcorr", action="store_true", dest="Zstack_pzIcorr",
                       default=0, help="intensity normalize piezo stacked images")
    parser.add_option("--fastpiezo", action="store_true", dest="fastpiezo",
                       default=0, help="piezo stack tiles are stored consecutively (instead of plane-wise)")
    #parser.add_option("--Ydown", action="store_true", dest="Ydown",
    #                   default=0, help="first pass is moving down (default is up)")
    parser.add_option("--short", action="store_const", const="short", dest="output_datatype",
                       default="byte", help="write short int data to file (default: byte)")
    parser.add_option("--scaleoutput",type="float",dest="scaleoutput", metavar="output scale", \
                       default=1.0, help="multiply slice images by scaleoutput before saving to file")
    parser.add_option("--file_type", type="string", dest="file_type", metavar="file_extension", \
                      default="tif",help="output file format (default: tif)")
    parser.add_option("--TV_file_type", type="string", dest="TV_file_type", metavar="file_extension", \
                      default="tif",help="TissueVision file format (default: tif)")
    parser.add_option("--use_IM", action="store_true", dest="im",
                       default=False, help="use imagemagick for preprocessing (old behaviour)")
    parser.add_option("--corr_tile_nonuniformity", action="store_true", dest="corr_tile_nonuniformity",
                       default=True, help="estimate and correct tile intensity nonuniformity")
    parser.add_option("--nocorr_tile_nonuniformity", action="store_false", dest="corr_tile_nonuniformity",
                       default=True, help="estimate and correct tile intensity nonuniformity")
    parser.add_option("--medfilter_tile", action="store_true", dest="medfilter_tile",
                       default=False, help="median filter the cropped tiles to eliminate 'spike' noise that cause spurious correlations")
    parser.add_option("--medfilter_size",type="int",dest="medfilter_size",default=3,\
                      help="size of median filter for cropped tiles to eliminate 'spike' noise")
    parser.add_option("--verbose", action="store_true", dest="verbose",
                       default=False, help="print output")
    parser.add_option("--keeptmp", action="store_true", dest="keeptmp",
                       default=False, help="keep temporary working directory (for debugging)")
    parser.add_option("--use_temp", type="string", dest="use_temp", metavar="file name", \
                      default=None,help="specify previous temp directory to use")
    parser.add_option("--skip_tile_match", action="store_true", dest="skip_tile_match",
                       default=False, help="skip tile matching and place tiles on perfect grid (for debugging)")

    options, args = parser.parse_args()
    VERBOSE = options.verbose

    try:
        if len(args)<2: raise FatalError("Specify input directory and output file.")
        inputdirectory = args[-2]
        outputfile = args[-1]
        if not options.clobber and os.path.exists(outputfile):
            raise FatalError("The --clobber option is needed to overwrite an existing file.")
    except FatalError as e:
        print('Error(%s):' % program_name, e.msg)
        raise SystemExit

    if (options.use_temp!=None):
        TEMPDIRECTORY=options.use_temp

    #generate image list
    starts=[]
    for j in [options.Zstart,options.Ystart,options.Xstart]:
        if (j>=0): starts.append(j)
        else: starts.append(None)
    ends=[]
    for j in [options.Zend,options.Yend,options.Xend]:
        if (j>=0): ends.append(j)
        else: ends.append(None)

    TileList,TVparamdict = generate_preprocessed_images(inputdirectory,starts=starts,ends=ends,\
                                                        channelflag=options.channel,imgftype=options.TV_file_type,\
                                                        fastpiezoloop=options.fastpiezo,gradcombine=options.gradimag,\
                                                        im=options.im,corr_tile_nonuniformity=options.corr_tile_nonuniformity,
                                                        medfilter_tile=options.medfilter_tile,medfilter_size=options.medfilter_size)
    uniqueZ=unique([ctile.indexarray[0] for ctile in TileList])

    #determine offsets with CCimages or read in positions from previously written file (or place images directly on a grid)
    existing_positions_file_flag = getattr(options,'use_positions_file')
    if (options.skip_tile_match):
        cmdout=run_subprocess("identify -format \"%%w %%h\" %s"%TileList[0].croppedfilename)
        matX=int(cmdout.split()[0]); matY=int(cmdout.split()[1])
        for ctile in TileList:
            ctile.pixoffsetarray[2] = ctile.indexarray[2]*matX
            ctile.pixoffsetarray[1] = ctile.indexarray[1]*matY
            ctile.pixoffsetarray[0] = ctile.indexarray[0]
    elif not existing_positions_file_flag:
        compute_offsets(TileList,overlapx=options.overlapx,overlapy=options.overlapy,Zref=options.Zref)
        if getattr(options,'save_positions_file'):
            save_positions_to_file(TileList,options.save_positions_file)
    else:
        get_positions_from_file(TileList,options.use_positions_file)

    #adjust positions to be all positive offsets based on global minima
    cmdout=run_subprocess("identify -format \"%%w %%h\" %s"%TileList[0].croppedfilename)
    matX=int(cmdout.split()[0]); matY=int(cmdout.split()[1])
    min_offset_x = min([ctile.pixoffsetarray[2] for ctile in TileList])
    min_offset_y = min([ctile.pixoffsetarray[1] for ctile in TileList])
    for ctile in TileList:
        ctile.pixoffsetarray[2] -= min_offset_x
        ctile.pixoffsetarray[1] -= min_offset_y
    max_offset_x = max([ctile.pixoffsetarray[2] for ctile in TileList])
    outimg_size_x = max_offset_x + matX
    max_offset_y = max([ctile.pixoffsetarray[1] for ctile in TileList])
    outimg_size_y = max_offset_y + matY

    #overlay images with opencv based run_image_overlay     
    Zstacklist=[]
    for z in uniqueZ:
        zinds = [cind for cind,ctile in enumerate(TileList) if ctile.indexarray[0]==z]
        clist = [TileList[j].croppedfilename for j in zinds]
        positions = array([TileList[j].pixoffsetarray for j in zinds],float)
        #generate full slices from tiles
        Zsliceimg=run_image_overlay(clist,positions[:,:],outimg_size_x=outimg_size_x,outimg_size_y=outimg_size_y,\
                                    outputfiletype=options.output_datatype,outscale=options.scaleoutput)
        Zstacklist.append(Zsliceimg)

    #if needed, perform slice-by-slice intensity normalization (i.e., for piezo stacks)
    if (options.Zstack_pzIcorr):
        j=0
        while (j*TVparamdict['N_z_piezo']<len(Zstacklist)):
            intensity_normalize_Zstack(Zstacklist[j*TVparamdict['N_z_piezo']:(j+1)*TVparamdict['N_z_piezo']],
                                       Zstacklist[j*TVparamdict['N_z_piezo']:(j+1)*TVparamdict['N_z_piezo']])
            j+=1

    #generate minc file
    if not (outputfile[-4:]==".mnc"): #output an image stack
        for z,cfile in enumerate(Zstacklist):
            cmdstr="cp %s %s"%(cfile,outputfile+'_Z%04d'%uniqueZ[z]+'.%s'%options.file_type)
            cmdout = run_subprocess(cmdstr)
    else: #output a mnc file
        x_step = 0.001*[TV_LORES,TV_HIRES][TVparamdict['mcolumns']>LORESMAT]
        y_step = 0.001*[TV_LORES,TV_HIRES][TVparamdict['mrows']>LORESMAT]
        if (TVparamdict['N_z_piezo']>1):              #this is wrong: can we incorporate both piezo and cut resolution in same mnc file??
            z_step = 0.001*2.0*TVparamdict['zres']    #2X for real resolution, apparently?!
            z_coord_list = 0.001*TVparamdict['sectionres']*(arange(len(Zstacklist))/TVparamdict['N_z_piezo']) + \
                           z_step*(arange(len(Zstacklist))%TVparamdict['N_z_piezo'])
            generate_mnc_file_from_tifstack(Zstacklist,outputfile,zstep=z_step,ystep=y_step,xstep=x_step,\
                                            outdatatype=options.output_datatype,Zcoordlist=z_coord_list) 
        else:
            z_step = 0.001*TVparamdict['sectionres']  
            generate_mnc_file_from_tifstack(Zstacklist,outputfile,zstep=z_step,ystep=y_step,xstep=x_step,outdatatype=options.output_datatype) 

    #clean up all temp files
    if (not options.keeptmp) and (options.use_temp==None):
        cmdout = run_subprocess("rm -r %s"%TEMPDIRECTORY)
    else:
        print("Temp directory is: %s"%TEMPDIRECTORY)



inputdirectory='/projects/souris/jgleave/TissueVision_2012/05Oct2012/Nestin_AF546'
TileList,TVparamdict = generate_preprocessed_images(inputdirectory,starts=[1,1,1],ends=[None,None,20])

from random import sample
from scipy.interpolate import bisplrep,bisplev
def avg_intensity_correction(TileList,ntiles=1024,scale=1.0):
    ##need to add some sort of selection process to avoid considering empty tiles
    #first average several tiles
    if (ntiles>len(TileList)):
        ilist = arange(len(TileList))
    else:
        ilist = sample(arange(len(TileList)),ntiles)
    cmdstr = "convert "
    for cindex in ilist:
        cmdstr += TileList[cindex].croppedfilename + ' '
    avgimage = gen_tempfile('iaverage',TileList[0].croppedfilename.split('.')[-1])
    cmdstr += "-level 0,%f%% -average %s"%(100.0/scale,avgimage)
    run_subprocess(cmdstr)
    #smooth average
    ifh = Image.open( avgimage )
    Iavg = reshape( asarray( ifh.getdata() ) ,(ifh.size[0],ifh.size[1]) )
    (Ny,Nx)=Iavg.shape
    stepsize=20
    Iavgdown=zeros(Iavg[stepsize-1::stepsize,stepsize-1::stepsize].shape,float)
    for j in range(stepsize):
        for k in range(stepsize):
            Iavgdown+=Iavg[j:stepsize*(Ny/stepsize):stepsize,k:stepsize*(Nx/stepsize):stepsize]
    Iavgdown /= float(stepsize**2)
    (Nyd,Nxd)=Iavgdown.shape
    tck=bisplrep(arange(Nyd*Nxd)%Nxd,arange(Nyd*Nxd)/Nxd,ravel(Iavgdown),task=0)
    z=transpose(bisplev(-0.5+arange(Ny)/float(stepsize),-0.5+arange(Nx)/float(stepsize),tck))
    z=z/mean(z.flat)
    Imul = where(z<0.5,2.0,1/z)
    #output Imul to file
    Imulfh = Image.fromarray(Imul)
    Imapimg = gen_tempfile('Imap',TileList[0].croppedfilename.split('.')[-1])
    Imulfh.save(Imapimg)
    #use to correct all tiles
    for ctile in TileList:
        tempname = ctile.croppedfilename[0:-4]+'_Iadj'+ctile.croppedfilename[-4:]
        cmdstr = 'convert %s %s -fx "u*v" %s'%(ctile.croppedfilename,Imapimg,tempname'
        run_subprocess(cmdstr)
        ctile.croppedfilename = tempname
    return 0

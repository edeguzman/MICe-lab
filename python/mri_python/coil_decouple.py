from sys import path as syspath
syspath.append('/home/bjnieman/source/mri_recon')
from mri_python.varian_read_file import *
from scipy.ndimage.morphology import binary_dilation
from scipy.linalg import lstsq,inv
from numpy.fft import *
from numpy.random import permutation
#from pylab import figure,plot,subplot,show,imshow,colorbar,cm

def fetch_grappa4prof_data(inputAcq,petable,remove_ppeshift=True,dcplppeadj=None): 
    #retrieve and reconstruct 4 2D "grappa profiles" from fid file
    nmice = inputAcq.nmice
    nacq = int(get_dict_value(inputAcq.param_dict,'np',1)/2)
    nro = int(get_dict_value(inputAcq.param_dict,'nro',1))
    nv = int(get_dict_value(inputAcq.param_dict,'nv',1))
    nv2 = int(get_dict_value(inputAcq.param_dict,'nv2',1))
    grappafov = int(get_dict_value(inputAcq.param_dict,'grappafov',1))
    nrcvrs = inputAcq.nrcvrs
    if (nmice>nrcvrs):
        nmice=nrcvrs
    t1_array = parse_petable_file(petable,'t1')
    t2_array = parse_petable_file(petable,'t2')
    mm_ppe = array(get_dict_value(inputAcq.param_dict,'mm_ppe',[0,0,0,0])).astype(float)
    lpe = get_dict_value(inputAcq.param_dict,'lpe',2.56)
    if (dcplppeadj==None):
        dcplppeadj = zeros((len(mm_ppe),),float)    
    #identify large fov data (even if not centred)
    pegrid = zeros((nv2*grappafov,nv*grappafov),int)
    pegrid[t2_array+(int(nv2/2)-1)*grappafov,t1_array+(int(nv/2)-1)*grappafov] = 1
    fovmask = zeros(pegrid.shape,int)
    fovmask[1:-1,1:-1] = ((pegrid[0:-2,0:-2] + pegrid[1:-1,0:-2] + pegrid[2::,0:-2] + \
                           pegrid[0:-2,1:-1] + pegrid[1:-1,1:-1] + pegrid[2::,1:-1] + \
                           pegrid[0:-2,2::]  + pegrid[1:-1,2::] + pegrid[2::,2::] ) > 2)*(pegrid[1:-1,1:-1]>0)
    i2,i1 = nonzero(fovmask)
    inds=array([],int)
    for j in range(len(i2)):
        inds = append(inds,nonzero( ((t2_array+(int(nv2/2)-1)*grappafov)==i2[j])&((t1_array+(int(nv/2)-1)*grappafov)==i1[j]) )[0] )
    grappapix = max(t2_array[inds])-min(t2_array[inds])+2
    grappapix = min(grappapix,max(t1_array[inds])-min(t1_array[inds])+2)                     
    raw_data = zeros((nmice,4,grappapix,nacq),complex)
    n_grabs = zeros((4,grappapix),int)
    t1cen = (max(t1_array[inds])+min(t1_array[inds]))/2; t1min = min(t1_array[inds])-int(mean(t1_array[inds]))
    t2cen = (max(t2_array[inds])+min(t2_array[inds]))/2; t2min = min(t2_array[inds])-int(mean(t2_array[inds]))
    if (remove_ppeshift):
        if (len(dcplppeadj)<len(mm_ppe)):
            dcplppeadj.extend([0 for x in range(len(mm_ppe)-len(dcplppeadj))])
        delppe = 2*pi*((mm_ppe+array(dcplppeadj,float))/lpe)
    else:
        delppe = zeros((nmice,),float)
    for imouse in range(nmice):
        for k in inds:
            fid_data,data_error = inputAcq.getdatafids(k,k+1,rcvrnum=imouse)
            if ( (t1_array[k]-t1cen)==0 ):
                raw_data[imouse,0,t2_array[k]-t2min,:]+=fid_data[0,:]
                if (imouse==0): n_grabs[0,t2_array[k]-t2min]+=1
            if ( (t2_array[k]-t2cen)==0 ):
                raw_data[imouse,1,t1_array[k]-t1min,:] += \
                                                   fid_data[0,:]*(exp(1.j*delppe[imouse]*float(t1_array[k])/float(grappafov))) 
                if (imouse==0): n_grabs[1,t1_array[k]-t1min]+=1
            if ( (t2_array[k]-t2cen)==(t1_array[k]-t1cen) ):
                raw_data[imouse,2,t2_array[k]-t2min,:] += \
                                                   fid_data[0,:]*(exp(1.j*delppe[imouse]*float(t1_array[k])/float(grappafov))) 
                if (imouse==0): n_grabs[2,t2_array[k]-t2min]+=1
            if ( (t2_array[k]-t2cen)==-(t1_array[k]-t1cen) ):
                raw_data[imouse,3,t2_array[k]-t2min,:] += \
                                                   fid_data[0,:]*(exp(1.j*delppe[imouse]*float(t1_array[k])/float(grappafov))) 
                if (imouse==0): n_grabs[3,t2_array[k]-t2min]+=1
    n_grabs = where(n_grabs>0,n_grabs,1)
    raw_data = raw_data/n_grabs[newaxis,:,:,newaxis]
    if (get_dict_value(inputAcq.param_dict,'sgflag','n')=='y'):
        maxind = raw_data.shape[-1]-nro+argmax(abs(raw_data[0,0,grappapix/2,-nro::]))
        endpt = [maxind+int(nro/2),nacq-1][maxind+int(nro/2)>=nacq]
        startpt = endpt-nro
    else:
        startpt = 0
        endpt = nro
        if (nro>nacq):
            print("ERROR: nro and np appear to be set inconsistently!!!")
            endpt = nacq
    profs = fftshift(fft2(fftshift(raw_data[:,:,:,startpt:endpt],axes=(-2,-1)),axes=(-2,-1)),axes=(-2,-1))
    return profs,0.5*(startpt+endpt),delppe


def mask_gappa4prof_data(maskshape,inputAcq):
    #profile directions/order established in fetch_grappa4prof_data
    #should really pass this or detect it
    prof_direc=array([[1,0,1],
                      [0,1,1],
                      [1,1,1],
                      [1,-1,1]],int)
    masks = zeros(maskshape,bool)
    #assume fov and mm_pro,mm_ppe and mm_ppe2 specify cylinder, orient specifies direction
    mm_ppe = array(get_dict_value(inputAcq.param_dict,'mm_ppe',[0,5.2,0,-5.4,-5.5,0,5.2,0,0,0,0,0,0,0,0,0])).astype(float)
    mm_ppe2 = array(get_dict_value(inputAcq.param_dict,'mm_ppe2',[0,-3,-6,-3,2.7,6,3,0,0,0,0,0,0,0,0,0,0])).astype(float)
    grappafov = get_dict_value(inputAcq.param_dict,'grappafov',8)
    lpe = get_dict_value(inputAcq.param_dict,'lpe',2.02)
    lpe2 = get_dict_value(inputAcq.param_dict,'lpe2',2.02)
    nmice = inputAcq.nmice
    nrcvrs = inputAcq.nrcvrs
    nmice = [nmice,nrcvrs][nmice>nrcvrs]
    ROIwidth = int( 1.5*maskshape[-2]/int(grappafov) )
    for j in range(nmice):
        for k in range(4): #4 profiles
            #fov and pixshift differences on diagonal axes taken care of by using unnormalized prof_direc in pixoffset calc
            pixoffset = maskshape[-2]*( (mm_ppe2[j]/(grappafov*lpe2))*prof_direc[k,0]+ 
                                        (mm_ppe[j] /(grappafov*lpe ))*prof_direc[k,1] )
            ROIinds = ( (pixoffset+maskshape[-2]/2+arange(ROIwidth)-ROIwidth/2)%maskshape[-2] ).astype(int)
            masks[j,k,ROIinds,:] = 1
    return masks


def find_axis_shift(imgdata1,maskdata1,imgdata2,maskdata2,start=-10,end=10,step=0.5,axis=-1):
    phaseramp=exp(1.j*(arange(imgdata1.shape[axis])-imgdata1.shape[axis]/2)*2*pi/imgdata1.shape[axis])
    steps=arange(start,end,step)
    cmask=(maskdata1!=maskdata2)
    if (not cmask.any()): cmask=ones(maskdata1.shape,bool)
    C = zeros((len(steps),),float)
    for j in range(len(steps)):
        pixshift = steps[j]
        imgdata1adj = fftshift(ifft(fftshift((phaseramp)**(pixshift/2.0)*
                      fftshift(fft(fftshift(imgdata1,axes=(axis,)),axis=axis),axes=(axis,)),axes=(axis,)),axis=axis),axes=(axis,))
        imgdata2adj = fftshift(ifft(fftshift((phaseramp)**(-pixshift/2.0)*
                      fftshift(fft(fftshift(imgdata2,axes=(axis,)),axis=axis),axes=(axis,)),axes=(axis,)),axis=axis),axes=(axis,))
        v1=abs(imgdata1adj)[cmask]
        v2=abs(imgdata2adj)[cmask]
        mv1=mean(v1)
        mv2=mean(v2)
        C[j] = sum((v1-mv1)*(v2-mv2))/sqrt(sum((v1-mv1)**2)*sum((v2-mv2)**2))
    i1 = argmax(C)
    istart = [0,i1-5][i1-5>0]
    iend = [len(steps),i1+5][i1+5<len(steps)]
    pfit=polyfit(steps[istart:iend],C[istart:iend],2)
    bestshift=-1*pfit[1]/(2.0*pfit[0])
    if ((bestshift>end) or (bestshift<start)):
        print("Read-out shift obtained is beyond search range...")
    return bestshift,polyval(pfit,bestshift)

def dcpl_axis_adjustment(profs,masks,cplgrps,nom_shift=10,axis=-1):
    #first, isolate read offsets (should be able to get this from mmrcvrf, get it here from files)
    ropos = zeros((profs.shape[0],),float)
    for cgrp in cplgrps:
        A = array([],float); b = array([],float);
        A = append(A,ones((len(cgrp),),float))
        b = append(b,0)
        for j in range(len(cgrp)-1):
            for k in range(j+1,len(cgrp)):
                shift,cval = find_axis_shift(profs[cgrp[j]],masks[cgrp[j]],profs[cgrp[k]],masks[cgrp[k]],
                                             start=-nom_shift,end=nom_shift,step=nom_shift/20.0,axis=axis)
                Arow = zeros((len(cgrp),),float)
                Arow[j] = cval; Arow[k] = -cval
                A = append(A,Arow)
                b = append(b,shift*cval)
        A.shape = (len(b),int(len(A)/len(b)))
        rofitpos,resids,rank,s = lstsq(A,b)
        ropos[cgrp] = rofitpos
    #apply offsets
    nro = profs.shape[axis]
    phaserampshape = ones(len(profs.shape)-1,int)
    phaserampshape[axis] = nro
    shifted_profs=empty(profs.shape,profs.dtype)
    for j in range(profs.shape[0]):
        phaseramp=exp(-1.j*(arange(nro)-nro/2)*2*pi*ropos[j]/nro)
        phaseramp.shape=phaserampshape
        kdata = fftshift(ifft(fftshift(profs[j,:,:,:],axes=(axis,)),axis=axis),axes=(axis,))
        shifted_profs[j,:,:,:] = fftshift(fft(fftshift( kdata*phaseramp ,axes=(axis,)),axis=axis),axes=(axis,))
    return shifted_profs,ropos
    
def find_coupling_coeffs(profs,masks,cplgrps,npts_lstsq=8000):
    #find coupling coeffs
    Cij = eye(profs.shape[0],dtype=complex)
    for cgrp in cplgrps:
        for j in range(len(cgrp)):
            cprofinds = list(set(cgrp)-set([cgrp[j]]))
            for k in cprofinds:
                cmask = (1-masks[cgrp[j]])*masks[k]*(sum(masks[cgrp],axis=0)<2)
                i3,i2,i1 = nonzero(cmask)
                if (len(i1)>npts_lstsq):
                    rinds = argsort(abs(profs[k])[i3,i2,i1])[-npts_lstsq::]
                    i3=i3[rinds]; i2=i2[rinds]; i1=i1[rinds]
                b=append(profs[cgrp[j]].real[i3,i2,i1],profs[cgrp[j]].imag[i3,i2,i1])
                A=zeros([2*len(i1),2],float)
                A[0:len(i1),0] = profs[k].real[i3,i2,i1]
                A[0:len(i1),1] = -profs[k].imag[i3,i2,i1]
                A[len(i1)::,0] = profs[k].imag[i3,i2,i1]
                A[len(i1)::,1] = profs[k].real[i3,i2,i1]
                cres,resid,rank,s = lstsq(A,b)
                Cij[cgrp[j],k] = cres[0]+1.j*cres[1]
    invCij = inv(Cij)
    print(invCij)
    return invCij

def output_grappa_profs(profs,dcpl_data,dcpl_profs_output_name,masks=None):
    profsmod = profs
    from pylab import figure,imshow,subplot,savefig
    for j in range(profs.shape[0]):
        figure(j+1);
        for k in range(4): #always 4 profiles
            if (masks!=None):
                edges = sqrt( abs(masks[j,k,:,:]-roll(masks[j,k,:,:],-1,axis=-2))**2+
                              abs(masks[j,k,:,:]-roll(masks[j,k,:,:],-1,axis=-1))**2 )
                maxsig = max(abs(profsmod[j,k,:,:]).flat)
                imgaspect = float(profs.shape[-1])/float(profs.shape[-2])
                if (imgaspect<0.6): imgaspect=0.6
                subplot(2,2,k+1); imshow(where(edges>0.5,maxsig,abs(profsmod[j,k,:,:])), 
                                         interpolation='nearest',aspect=imgaspect)
            else:
                subplot(2,2,k+1); imshow(abs(profsmod[j,k,:,:]))
        c_outname = dcpl_profs_output_name[:-4]+"_%d"%j+dcpl_profs_output_name[-4:]
        savefig(c_outname)
        print("Output decoupling profile for coil#%d to %s..."%(j,c_outname))
    return None

class decouple_info:
    def __init__(self,roposition,invCij,rok0index,nro,ppeshift):
        self.roposition=roposition
        self.invCij=invCij
        self.rok0index=int(rok0index)
        self.nro=int(nro)
        self.ppeshift=ppeshift

def gen_decoupling_info(inputAcq,petable,cplgrps_string="0,2,5;1,3,4,6",
                        dcplppeadj=None,remove_ppeshift=True,dcpl_profs_output_name=None):
    print('Generating coil decoupling coefficients and offsets...')
    profs,rok0index,delppe = fetch_grappa4prof_data(inputAcq,petable,remove_ppeshift=remove_ppeshift,dcplppeadj=dcplppeadj)
    masks = mask_gappa4prof_data(profs.shape,inputAcq)
    cplgrps = [[int(x) for x in grpstrings.split(',')] for grpstrings in cplgrps_string.split(';')]
    if (max([max(x) for x in cplgrps])>=profs.shape[0]):
        print("Inconsistency between grappa_coil_groupings and acquired channels...")
        raise SystemExit
    profs,ropos = dcpl_axis_adjustment(profs,masks,cplgrps,nom_shift=int( inputAcq.param_dict['sw']/2000.0 ),axis=-1) 
    invCij = find_coupling_coeffs(profs,masks,cplgrps)
    dcpl_data = decouple_info(ropos,invCij,rok0index,profs.shape[-1],delppe) 
    if (dcpl_profs_output_name!=None):
        output_grappa_profs(profs,dcpl_data,dcpl_profs_output_name,masks=masks)    
    return dcpl_data





import re
import numpy as N
from numpy.fft import fftshift,fft
from numpy.fft import fft2 as fft2d
from numpy.fft import fftn as fftnd

#from LinearAlgebra import *
import mri_python.varian_read_file as vrf

from scipy.signal import medfilt

########################################################################################################################################################    
#simple function to fetch parameter with a default if absent in procpar file

class FatalError(Exception):
    def __init__(self,args=None):
        self.msg = args

########################################################################################################################################################    
#simple function to fetch parameter with a default if absent in procpar file

def get_dict_value(image_param_dict,key,default):
    retvalue=default
    if key in image_param_dict:
        retvalue=image_param_dict[key]
    return retvalue

########################################################################################################################################################    
#k-space reordering functions

def centre_out_reordering(data,axis=-2):
    print("Reordering centre-out acquisition order (axis %d)..." % axis)
    dim_sizes=N.shape(data)
    indices=range(dim_sizes[axis]-1,0,-2)+range(0,dim_sizes[axis],2)
    reorder_data=N.take(data,indices,axis)
    return reorder_data

def petable_reordering(data,axis=-2,petable_array='t1',petable_name=''):
    print("Reordering from petable file (axis %d)..." % axis)
    dim_sizes=N.shape(data)
    enc_array = vrf.parse_petable_file(petable_name,petable_array)
    enc_array = enc_array - N.minimum.reduce(enc_array) #this is used to avoid 0-N-1 vs. 1-N differences
    indices = N.zeros((dim_sizes[axis],),N.int)
    N.put(indices,enc_array,N.arange(dim_sizes[axis]))
    reorder_data=N.take(data,indices,axis)
    return reorder_data

def petable_orderedpair_reordering(data,petable_arrays=('t1','t2'),petable_name='',matrix=(None,None),index_start=0,index_end=None,t1array=None,t2array=None):
    dimnums=(-2,-3) #only handles this arrangement for now, any reason to do otherwise?
    print("Ordered pair reordering from petable file (axes %d, %d)..." % (dimnums[0],dimnums[1]))
    dim_sizes=N.shape(data)
    if ((t1array is None) and (t2array is None)):
        enc_array_1 = vrf.parse_petable_file(petable_name,petable_arrays[0])
        enc_array_1 = enc_array_1 - N.minimum.reduce(enc_array_1)
        enc_array_2 = vrf.parse_petable_file(petable_name,petable_arrays[1])
        enc_array_2 = enc_array_2 - N.minimum.reduce(enc_array_2)
    else:
        enc_array_1 = t1array - N.minimum.reduce(t1array)
        enc_array_2 = t2array - N.minimum.reduce(t2array)
    nv = [matrix[0],1+N.maximum.reduce(enc_array_1)][matrix[0]==None]
    nv2 = [matrix[1],1+N.maximum.reduce(enc_array_2)][matrix[1]==None]    
    outputdata_shape = N.array(data.shape)
    outputdata_shape[-3] = nv2
    outputdata_shape[-2] = nv
    datashape = N.array(data.shape); datashape[-3]*=datashape[-2]; datashape=N.delete(datashape,-2)
    data = N.reshape(data,tuple(datashape))
    outdata = N.zeros( outputdata_shape , data.dtype)
    nacq = N.zeros( (nv2,nv),int)
    if (not index_end):
        index_end = len(enc_array_1)
    enc_array_1 = enc_array_1[index_start:index_end]
    enc_array_2 = enc_array_2[index_start:index_end]
    #outdata[...,enc_array_2,enc_array_1,:] = data
    for j in range(data.shape[-2]):
        j_encind = j%len(enc_array_1)
        outdata[...,enc_array_2[j_encind],enc_array_1[j_encind],:] += data[...,j,:]
        nacq[enc_array_2[j_encind],enc_array_1[j_encind]]+=1
    nacq = N.where(nacq<1,1,nacq)
    return outdata/nacq[:,:,N.newaxis]

########################################################################################################################################################    
#FOV shifting functions

def fov_shift(data,pixel_shift,dimnum,large_data_flag):
    dim_sizes = data.shape
    print("Performing FOV shift (axis %d, %f of %d pixels)..." % (dimnum,pixel_shift,dim_sizes[dimnum]))
    swapped_data=N.swapaxes(data,dimnum,-1)
    if (type(pixel_shift)!=float) and (type(pixel_shift)!=N.float64): pixel_shift=pixel_shift[0]
    if (large_data_flag):
        swap_shape = swapped_data.shape
        phase_ramp= (N.exp(-N.arange(dim_sizes[dimnum],dtype=N.float16)*comp_i*2*N.pi*pixel_shift/dim_sizes[dimnum]))
        phase_ramp *= N.conjugate(phase_ramp[dim_sizes[dimnum]/2])
        for j in range(swap_shape[0]):
            swapped_data[j] = (swapped_data[j]*phase_ramp).astype(N.complex)
    else:
        phase_ramp=N.exp(-N.arange(dim_sizes[dimnum],dtype=N.float16)*1.j*2*N.pi*pixel_shift/dim_sizes[dimnum])
        phase_ramp *= N.conjugate(phase_ramp[dim_sizes[dimnum]//2])
        swapped_data=swapped_data*phase_ramp
    shifted_data=N.swapaxes(swapped_data,dimnum,-1)
    return shifted_data

def fov_adjustment(rawdata,options,inputAcq,mouse_num):
    #ro shift, if specified on command line
    if hasattr(options,'fov_shift_ro'):
        ro_pixshift = float(options.fov_shift_ro.split(',')[mouse_num])
        if (abs(ro_pixshift)>0.1):
            print("RO shift...")
            rawdata=fov_shift(rawdata,ro_pixshift,-1,options.large_data_recon)
    #pe1 shift, by file and by user specification
    if (inputAcq.platform=="Varian"):
        prescribed_pe1_shift = get_dict_value(inputAcq.param_dict,'mm_ppe',[0.0,0.0,0.0,0.0])[mouse_num]
        pe1_pix = float( get_dict_value(inputAcq.param_dict,'lpe',14.0) )/get_dict_value(inputAcq.param_dict,'nv',432)
    elif (inputAcq.platform=="Bruker"):
        #need to incorporate offsets from inputAcq.SliceGeo['Offsets'] and from MICe_*_offset 
        prescribed_pe1_shift = get_dict_value(inputAcq.method_param_dict,'MICe_Phase1Offset',[0.0,0.0,0.0,0.0])[mouse_num]
        pe1_pix = float( get_dict_value(inputAcq.method_param_dict,'PVM_Fov',[14.0,14.0,14.0])[1] )/inputAcq.npe   
    pe1_pixshift = prescribed_pe1_shift/pe1_pix
    if (options.noshift_ppe): pe1_pixshift=0.0
    if hasattr(options,'fov_shift_pe1'):
        pe1_pixshift += float(options.fov_shift_pe1.split(',')[mouse_num])
    if (abs(pe1_pixshift)>0.01):
        rawdata=fov_shift(rawdata,pe1_pixshift,-2,options.large_data_recon)
    #pe2 shift if required, by file and by user specification    
    if (inputAcq.nD>2): #check for 2D vs. 3D should be here instead of default run
        if (inputAcq.platform=="Varian"):
            prescribed_pe2_shift = get_dict_value(inputAcq.param_dict,'mm_ppe2',[0.0,0.0,0.0,0.0])[mouse_num]
            pe2_pix = float( get_dict_value(inputAcq.param_dict,'lpe2',14.0) )/get_dict_value(inputAcq.param_dict,'nv2',432)
        elif (inputAcq.platform=="Bruker"):
            #prescribed_pe2_shift = get_dict_value(inputAcq.method_param_dict,'MICe_Phase2Offset',[0.0,0.0,0.0,0.0])[mouse_num]
            prescribed_pe2_shift = get_dict_value(inputAcq.method_param_dict,'MICe_SliceOffset',[0.0,0.0,0.0,0.0])[mouse_num]
            pe2_pix = float( get_dict_value(inputAcq.method_param_dict,'PVM_Fov',[14.0,14.0,14.0])[2] )/inputAcq.npe2   
        pe2_pixshift = prescribed_pe2_shift/pe2_pix
        if (options.noshift_ppe2): pe2_pixshift=0.0
        if hasattr(options,'fov_shift_pe2'):
            pe2_pixshift += float(options.fov_shift_pe2.split(',')[mouse_num])
        if (abs(pe2_pixshift)>0.5):
            rawdata=fov_shift(rawdata,pe2_pixshift,-3,options.large_data_recon)
    return rawdata

########################################################################################################################################################    
#simplest version of fetching k-space data from fid file

def gen_kspace_simple(inputAcq,ircvr):
    print("Generating simple k-space data...")
    nrcvrs = inputAcq.nrcvrs
    no_ro = inputAcq.nro
    if ((no_ro<=0) or (no_ro>inputAcq.data_shape[-1])): no_ro = inputAcq.data_shape[-1] #np
    if (len(inputAcq.data_shape)>2) and (inputAcq.data_shape[-3]>0):
        npe2 = inputAcq.data_shape[-3]
    else:
        npe2 = 1
    npe1 = inputAcq.data_shape[-2]
    if (inputAcq.platform=="Varian"):
        nf = inputAcq.nf
        ni_perchan = inputAcq.nfid*inputAcq.ni
        nreps = 1 
    elif (inputAcq.platform=="Bruker"): #force Bruker to artificial Varian format
        nf = inputAcq.nf #inputAcq.npe
        nreps = get_dict_value(inputAcq.method_param_dict,"PVM_NRepetitions",1)
        if (inputAcq.nD==2):
            ni_perchan = inputAcq.nslices
        elif (inputAcq.nD==3):        
            ni_perchan = inputAcq.ni_perchan #npe2
    else:
        print("Input format not recognized.")
    raw_imgdata = N.zeros((nreps,ni_perchan,nf,no_ro),N.complex)
    for k in range(nreps):
        for j in range(ni_perchan):
            print("%d / %d" % (k*ni_perchan+j,nreps*ni_perchan))
            fid_start = j*nf
            fid_end = (j+1)*nf
            fid_data,data_error = inputAcq.getdatafids(fid_start,fid_end,rcvrnum=ircvr)
            if (data_error):
                data_fraction = float(j)/float(ni_perchan)
                print('Error at %f%% through data set' % (100.0*data_fraction))
                break
            raw_imgdata[k,j,:,:] = fid_data[:,-no_ro:]
    if ((inputAcq.nD==2) and (inputAcq.nslices>1)): #for interleaved slices
        raw_imgdata_sorted = N.zeros(raw_imgdata.shape,raw_imgdata.dtype)
        raw_imgdata.shape = (raw_imgdata.shape[-3]*raw_imgdata.shape[-2],raw_imgdata.shape[-1])
        print(raw_imgdata.shape)
        for j in range(raw_imgdata_sorted.shape[-3]):
            raw_imgdata_sorted[...,j,:,:] = raw_imgdata[j::inputAcq.nslices,:]
        raw_imgdata = raw_imgdata_sorted
    if (len(inputAcq.data_shape)==5) and (inputAcq.data_shape[0]>0):
        try: 
            raw_imgdata = N.reshape(raw_imgdata,(inputAcq.data_shape[0],ni_perchan,nf,no_ro))
        except ValueError:
            raw_imgdata = N.reshape(raw_imgdata,(ni_perchan,nf,no_ro))
    else:
        raw_imgdata = N.reshape(raw_imgdata,(ni_perchan,nf,no_ro))
    return raw_imgdata

def DCartcorr(kdata,param_dict):
    (npe2,npe,nro)=kdata.shape
    img=fftshift(fftnd(fftshift(kdata)))
    maxind = N.argmax(N.ravel(abs(img)))
    (pe2ind,pe1ind,roind) = (maxind/(npe*nro),(maxind/nro)%npe,maxind%nro)
    direc = [1,-1][int(abs(img[pe2ind,pe1ind-1,roind])>abs(img[pe2ind,pe1ind+1,roind]))]
    delta=direc * abs(img[pe2ind,pe1ind+direc,roind]) / ( abs(img[pe2ind,pe1ind,roind]) + abs(img[pe2ind,pe1ind+direc,roind]) )
    pos = pe1ind+delta-npe/2
    direc = [1,-1][int(abs(img[pe2ind,pe1ind,roind-1])>abs(img[pe2ind,pe1ind,roind+1]))]
    delta=direc * abs(img[pe2ind,pe1ind,roind+direc]) / ( abs(img[pe2ind,pe1ind,roind]) + abs(img[pe2ind,pe1ind,roind+direc]) )
    posro = roind+delta-nro/2
    kshift = kdata* (N.exp(-1.j*2*N.pi*pos*(N.arange(npe)-npe/2)/float(npe))[N.newaxis,:,N.newaxis]) * \
                    (N.exp(-1.j*2*N.pi*posro*(N.arange(nro)-nro/2)/float(nro))[N.newaxis,N.newaxis,:])
    imgshift = fftshift(fftnd(fftshift(kshift)))
    dcoff = (imgshift[npe2/2,npe/2,nro/2]-0.5*(imgshift[npe2/2,npe/2-1,nro/2]+imgshift[npe2/2,npe/2+1,nro/2]))/(nro*npe*npe2)
    kdatacorr = -1.0*dcoff*N.exp(1.j*2*N.pi*pos*(N.arange(npe)-npe/2)/float(npe))[:,N.newaxis]*N.exp(1.j*2*N.pi*posro*(N.arange(nro)-nro/2)/float(nro))[N.newaxis,:]
    kdatacorr = kdatacorr[N.newaxis,:,:] + kdata
    return kdatacorr,dcoff

########################################################################################################################################################    
#k-space filtering functions

def fermi_ellipse_filter(rawdata,cutoff=0.98,transwidth=0.06):
    print("Applying elliptical fermi filter...")
    data_shape = N.shape(rawdata)
    nro = data_shape[-1]
    nv1 = data_shape[-2]
    nv2 = data_shape[-3]
    for j in range(nv2):
        r_vals = N.sqrt( (2.0*j/float(nv2)-1.0)**2.0 + \
                       ((2.0*N.arange(nv1)/float(nv1)-1.0)**2.0)[:,N.newaxis] + \
                       ((2.0*N.arange(nro)/float(nro)-1.0)**2.0)[N.newaxis,:]  \
                     )
        filt = 1.0/(1.0+N.exp(-(cutoff-r_vals)/transwidth))
        if len(data_shape)==3:
            rawdata[j,:,:]=(rawdata[j,:,:]*filt).astype(N.complex)
        else:
            rawdata[...,j,:,:]=(rawdata[...,j,:,:]*filt).astype(N.complex)
        #if len(data_shape)==4:
        #    for k in range(data_shape[0]):
        #        rawdata[k,j,:,:] = (rawdata[k,j,:,:]*filt).astype(N.complex)
        #elif len(data_shape)==3:
        #    rawdata[j,:,:] = (rawdata[j,:,:]*filt).astype(N.complex)
    return rawdata

def fermi_pecircle_filter(rawdata):
    print("Applying circular pe fermi filter...")
    cutoff = 0.98
    transwidth = 0.06
    data_shape = N.shape(rawdata)
    nv1 = data_shape[-2]
    nv2 = data_shape[-3]
    r_vals = N.sqrt( ((2.0*N.arange(nv2)/float(nv2)-1.0)**2.0)[:,N.newaxis] + \
                   ((2.0*N.arange(nv1)/float(nv1)-1.0)**2.0)[N.newaxis,:] \
                 )
    filt = 1.0/(1.0+N.exp(-(cutoff-r_vals)/transwidth))
    if len(data_shape)==4:
        for k in range(data_shape[0]):
            rawdata[k,:,:,:] = (rawdata[k,:,:,:]*filt[:,:,N.newaxis]).astype(N.complex)
    elif len(data_shape)==3:
        rawdata[:,:,:] = (rawdata[:,:,:]*filt[:,:,N.newaxis]).astype(N.complex)
    return rawdata


########################################################################################################################################################    
#FFT functions for recon
def recon_2d_slices(rawdata):
    print("Performing 2D FFT reconstruction...")
    image_data = fftshift( fft2d( fftshift(rawdata,axes=(-2,-1)) ,s=None,axes=(-2,-1)),axes=(-2,-1))
    return image_data

def recon_3d(rawdata):
    print("Performing 3D FFT reconstruction...")
    image_data = fftshift( fftnd( fftshift(rawdata,axes=(-3,-2,-1)) ,s=None,axes=(-3,-2,-1)),axes=(-3,-2,-1))
    return image_data

def recon_1D_RO(rawdata):
    print("Performing 1D FFT...")
    image_data = fftshift( fft( fftshift(rawdata,axes=(-1,)) ,axis=-1),axes=(-1,))
    return image_data

########################################################################################################################################################    
#default handling of FFT

def default_recon(seqrec):
    if seqrec.options.nofft:
        print("No fft...")
        seqrec.image_data=seqrec.kspace
    elif seqrec.options.fft2d:
        seqrec.image_data=recon_2d_slices(seqrec.kspace)
    elif seqrec.options.fft1d:
        seqrec.image_data=recon_1D_RO(seqrec.kspace)
    else:
        sizes = seqrec.kspace.shape
        if (len(sizes)==4):
            seqrec.image_data=N.zeros(sizes,N.complex)
            for j in range(sizes[0]):
                seqrec.image_data[j,:,:,:] = recon_3d(seqrec.kspace[j,:,:,:]).astype(N.complex)
        elif (seqrec.options.large_data_recon):
            print("Large data (in place) reconstruction...")
            if (len(sizes)==3):
                for j in range(sizes[-1]):
                    seqrec.kspace[:,:,j] = fftshift( fft2d( fftshift(seqrec.kspace[:,:,j]
                                                           ,axes=(0,1)) ,s=None,axes=(-2,-1)),axes=(0,1)).astype(N.complex)
                for j in range(sizes[0]):
                    seqrec.kspace[j,:,:] = fftshift( fft( fftshift(seqrec.kspace[j,:,:],axes=(-1,)),axis=-1),axes=(-1,)).astype(N.complex)
                seqrec.image_data = seqrec.kspace
        else:
            seqrec.image_data=recon_3d(seqrec.kspace).astype(N.complex)

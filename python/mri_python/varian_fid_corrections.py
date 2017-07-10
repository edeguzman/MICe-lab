from scipy.interpolate import LSQUnivariateSpline
from scipy.stats import mode
import collections
from mri_python.varian_read_file import *
from itertools import combinations
from numpy.fft import *
from scipy.linalg import lstsq

def phase_drift_corr(inputAcq,petable,imouse=None,petable_arrays=('t1','t2')):
    if (inputAcq.platform!="Varian"):
        "Function phase_drift_corr only functions for Varian acquisitions..."
        raise SystemExit
    print('Estimating smoothed phase drift correction...')
    #get repeated k0 grabs
    nacq = int(get_dict_value(inputAcq.param_dict,'np',1)/2)
    nro = int(get_dict_value(inputAcq.param_dict,'nro',1))
    etl = int(get_dict_value(inputAcq.param_dict,'etl',1))
    nrcvrs = inputAcq.nrcvrs
    nmice = inputAcq.nmice
    if (nmice>nrcvrs):
        nmice=nrcvrs
    t1_array = parse_petable_file(petable,petable_arrays[0])
    t2_array = parse_petable_file(petable,petable_arrays[1])
    i1 = nonzero( (t1_array==0)*(t2_array==0) )[0]
    if (len(i1)<20):
        print('Too few k0 grabs for phase drift correction...')
        return ones(len(t1_array),float)
    if (etl>1):
        print('Phase drift correction not ready for etl>1...')
        return zeros(len(t1_array),float)
    if (imouse==None):
        mouselist = range(nmice)
    else:
        mouselist = [imouse]
    k0_data = zeros((len(mouselist), len(i1), nacq), complex)
    for k in range(len(mouselist)):
        for j in range(len(i1)):
            fid_data,data_error = inputAcq.getdatafids(i1[j],i1[j]+1,rcvrnum=k)
            k0_data[k,j,:]=fid_data.copy()
    #evaluate phase drift from k0 grabs
    maxind = k0_data.shape[-1]-nro+argmax(abs(k0_data[0,0,-nro:]),axis=-1)
    phasecorr = zeros((len(mouselist),len(t2_array)),complex)
    Nacqs = len(t2_array)
    i1mod = append(i1,len(t2_array)) 
    for j in range(k0_data.shape[0]):
        refangles = unwrap(angle( k0_data[j,:,maxind] ))
        refangles = append(refangles,median(refangles[-11::])) # artificially add a last point to limit 
                                                               # behaviour at the end
        splsmooth = LSQUnivariateSpline(i1mod,refangles,t=arange(Nacqs/3,int(0.8*Nacqs),Nacqs/3))
        smoothphase = splsmooth(arange(Nacqs))
        smoothphase = median(smoothphase)-smoothphase #this flips the sign
        phasecorr[j,:] = exp(1.j*smoothphase)
    return phasecorr


def get_corrected_datafids(inputAcq,fid_start,fid_end,mouse_num=0,phasecorr=None,dcpl_info=None,dcpl_ppe_index=0):
    if (inputAcq.platform!="Varian"):
        "Function get_corrected_datafids only functions for Varian acquisitions..."
        raise SystemExit
    if ((dcpl_info is None) and (phasecorr is None)): #no corrections to apply
         fid_data,data_error = inputAcq.getdatafids(fid_start,fid_end,rcvrnum=mouse_num)
    elif (dcpl_info is None) and (phasecorr is not None):
        fid_data,data_error = inputAcq.getdatafids(fid_start,fid_end,rcvrnum=mouse_num)
        fid_data = fid_data*phasecorr[mouse_num,fid_start:fid_end,newaxis]
    else:
        cgrp = nonzero( abs(dcpl_info.invCij[mouse_num,:])>1e-4 )[0]
        #start = dcpl_info.rok0index-dcpl_info.nro/2
        #end = start+dcpl_info.nro
        np = inputAcq.header_info[2]
        fid_data = zeros((fid_end-fid_start,int(np/2)),complex)  #dcpl_info.nro
        for j in cgrp:
            cfid,data_error = inputAcq.getdatafids(fid_start,fid_end,rcvrnum=j)
            fid_data += cfid[:,:]* \
                        dcpl_info.invCij[mouse_num,j]* \
                        exp(-1.j*2*pi*(arange(cfid.shape[-1])-dcpl_info.rok0index)*dcpl_info.roposition[j]/dcpl_info.nro)* \
                        exp(1.j*dcpl_info.ppeshift[j]*dcpl_ppe_index)
        if (phasecorr is not None):
            fid_data = fid_data*phasecorr[mouse_num,fid_start:fid_end,newaxis] 
    if (data_error):
        print('Unable to retrieve all fids (%d,%d)...' % (fid_start,fid_end))
    return fid_data



def rep_pos_corr(inputAcq,petable,imouse=None,petable_arrays=('t1','t2'),corrmat=64,phasedriftcorr=None):
    if (inputAcq.platform!="Varian"):
        "Function get_corrected_datafids only functions for Varian acquisitions..."
        raise SystemExit
    print('Estimating position shift between reps...')
    nacq = int(get_dict_value(inputAcq.param_dict,'np',1))/2
    nro = int(get_dict_value(inputAcq.param_dict,'nro',1))
    etl = int(get_dict_value(inputAcq.param_dict,'etl',1))
    nv = int(get_dict_value(inputAcq.param_dict,'nv',1))
    nv2 = int(get_dict_value(inputAcq.param_dict,'nv2',1))
    grappafov = int(get_dict_value(inputAcq.param_dict,'grappafov',1))
    t1_array = parse_petable_file(petable,petable_arrays[0])
    t2_array = parse_petable_file(petable,petable_arrays[1])
    i1 = nonzero( (t1_array/grappafov>-corrmat/2)*(t1_array/grappafov<=corrmat/2)* \
                  (t2_array/grappafov>-corrmat/2)*(t2_array/grappafov<=corrmat/2)* \
                  (t1_array%grappafov==0)*(t2_array%grappafov==0) )[0]
    noutreps = int( mode(array( collections.Counter(t1_array[i1]+nv/2-1+nv*(t2_array[i1]+nv2/2-1)).most_common() )[:,1])[0][0] )
    if (etl>1):
        print('Correction not ready for etl>1...')
        return ones(len(t1_array),float)
    if (imouse==None):
        mouselist = range(inputAcq.nmice)
    else:
        mouselist = [imouse]
    if (phasedriftcorr==None):
        phasedriftcorr=ones((len(mouselist),len(t1_array)),float)
    kdata = zeros((len(mouselist),noutreps,corrmat,corrmat,nacq),complex)
    for k in range(len(mouselist)):
        for q in range(corrmat):
            for r in range(corrmat):
                inds = i1[ nonzero( (t1_array[i1]==(r-corrmat/2+1)*grappafov)*(t2_array[i1]==(q-corrmat/2+1)*grappafov) )[0] ]
                repstep = len(inds)/noutreps
                for j in range(len(inds)):
                     fid_data,data_error = inputAcq.getdatafids(fid_start,fid_end,rcvrnum=k)
                     kdata[k,[noutreps-1,j/repstep][j/repstep<noutreps],q,r,:]+=fid_data[0,:]*phasedriftcorr[k,inds[j]]/repstep
    #evaluate phase ramp drift
    kplanes = zeros((len(mouselist),noutreps,corrmat,corrmat),complex)
    for k in range(len(mouselist)):
        if (get_dict_value(inputAcq.param_dict,'sgflag','n')=='y'):
            maxind = kdata.shape[-1]-nro+argmax(abs(kdata[k,0,corrmat/2,corrmat/2,-nro:]),axis=-1)
        else:
            maxind = argmax(abs(kdata[k,0,corrmat/2,corrmat/2,:]),axis=-1)
        kplanes[k,:,:,:] = kdata[k,:,:,:,maxind]
    repstep = len(t1_array)/noutreps
    repind = arange(len(t1_array))/repstep
    repind = where(repind>=noutreps,noutreps-1,repind)
    for j in range(len(mouselist)):
        p=zeros((2*noutreps,),float)
        peshifts=PEadjustment(kplanes[j,:,:,:])
        delphase = exp(1.j*2*pi*(arange(nv2)-nv2/2)[newaxis,:]*peshifts[:,0,newaxis]/corrmat)[:,:,newaxis]* \
                   exp(1.j*2*pi*(arange(nv)-nv/2)[newaxis,:]*peshifts[:,1,newaxis]/corrmat)[:,newaxis,:]
        phasedriftcorr[j,:] = phasedriftcorr[j,:]*delphase[repind,t2_array/grappafov+nv2/2-1,t1_array/grappafov+nv/2-1]
        print("PE2 pixel shifts (image %d): "%mouselist[j],peshifts[:,0]*float(nv2)/float(corrmat))
        print("PE1 pixel shifts (image %d): "%mouselist[j],peshifts[:,1]*float(nv)/float(corrmat))
    return phasedriftcorr    


def PEshift(kplane1,kplane2,retpixshift=True):
    # need to select planes instead of lines
    angdiff = angle(kplane1*conj(kplane2))
    A=ones((prod(angdiff.shape),3),float)
    A[:,0]=arange(prod(angdiff.shape))/angdiff.shape[1]-angdiff.shape[0]/2
    A[:,1]=arange(prod(angdiff.shape))%angdiff.shape[1]-angdiff.shape[1]/2    
    x,resids,rank,s = lstsq(A,ravel(angdiff))
    if (retpixshift):
        return x[0]*angdiff.shape[0]/(2*pi),x[1]*angdiff.shape[1]/(2*pi)
    else:
        return x


def PEadjustment(kmdata):
    pepos = zeros((kmdata.shape[0],2),float)
    A = array([],float); b = array([],float); w = array([],float)
    for j in range(2):
        Arow = zeros((2*kmdata.shape[0],),float); Arow[j::2] = 1.0
        A = append(A,Arow)
        b = append(b,0)
    for q in combinations(arange(kmdata.shape[0]),2):
        pe2shift,pe1shift = PEshift(kmdata[q[0]],kmdata[q[1]])
        Arow = zeros((2*kmdata.shape[0],),float); 
        Arow[2*q[0]] = 1; Arow[2*q[1]] = -1
        A = append(A,Arow); b = append(b,pe2shift); 
        Arow = zeros((2*kmdata.shape[0],),float); 
        Arow[2*q[0]+1] = 1; Arow[2*q[1]+1] = -1
        A = append(A,Arow); b = append(b,pe1shift); 
    A.shape = (len(b),len(A)/len(b))
    pefitpos,resids,rank,s = lstsq(A,b)
    pefitpos.shape=(kmdata.shape[0],2)
    return -1*pefitpos

    

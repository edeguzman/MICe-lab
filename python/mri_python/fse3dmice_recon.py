#!/usr/bin/env python
#program to process the fse3dmice_lsn_gshape_cyl_loop MRI experiment.
#August, 2014
#Brian Nieman and Leigh Spencer Noakes


import os
from sys import path as syspath
from optparse import OptionGroup
syspath.append('/home/bjnieman/source/mri_recon') 
from mri_python.recon_genfunctions import *  
import mri_python.varian_read_file as vrf
#from varian_recon import *
from pylab import figure,plot,subplot,show,imshow,colorbar,cm
from scipy.optimize import leastsq
from numpy.linalg import lstsq
from numpy import *
from numpy.fft import*



def seq_specific_options(parser):
    optgroup = OptionGroup(parser,"CYLfse3dmice","sequence-specific reconstruction options")
    optgroup.add_option("--phasecorr_data",type="string",dest="phasecorr_data",
                       default=None, help="path to correction fid data")
    optgroup.add_option("--phasecorr_table",type="string",dest="phasecorr_table",
                       default=None, help="petable for phasecorr")
    optgroup.add_option("--no_even_odd_echo_roshift", action="store_true", dest="no_even_odd_echo_roshift",
                        default=0, help="DOES NOT perform readout k-pos shift between even and odd echoes based on phasecorr data")
    optgroup.add_option("--no_even_odd_echo_phaseshift", action="store_true", dest="no_even_odd_echo_phaseshift",
                        default=0, help="DOES NOT perform readout phase shift between even and odd echoes based on phasecorr data")
    optgroup.add_option("--no_moment_terms_from_petable", action="store_true", dest="no_moment_terms_from_petable",
                        default=0, help="DOES NOT calculates moments for pe1 pe2 pe1^2 pe2^2 pe1pe2 (x,z,xz phase corrections)")
    optgroup.add_option("--no_Echo_shift_apply", action="store_true", dest="no_Echo_shift_apply",
                        default=0, help="DOES NOT apply echo shift based on centre echo")
    optgroup.add_option("--echoamp_alpha",type="float",dest="echoamp_alpha",
                      default=0.02, help="alpha in echo amplitude correction term (default 0.02): (Ei/E0)/( (Ei/E0)**2 + alpha )") 
    parser.add_option_group(optgroup)              
  
               
    
class seq_reconstruction():
    def __init__(self,inputAcq,options,outputfile,noshift_ppe=True,noshift_ppe2=False):
        if (inputAcq.platform!="Varian"):
            print("fse3dmice_recon is for Varian acquisitions only!!")
            raise SystemExit
        self.options = options
        self.inputAcq = inputAcq
        self.outputfile = outputfile
        if (self.options.petable==''):
            self.petable_name = vrf.get_dict_value(inputAcq.param_dict,'petable','petable')
        else:
            self.petable_name = self.options.petable
        self.kspace=None
        #create all arrays
        self.alpha = self.options.echoamp_alpha #coefficient for amplitude correction filter
        self.even_odd_pixel_shift = zeros(inputAcq.nmice,float)
        self.ephase_corr = ones((inputAcq.nmice,vrf.get_dict_value(inputAcq.param_dict,'etl',6)),complex)
        self.pMfit = zeros((inputAcq.nmice,5),float)
        self.maxroind=None 
        self.options.noshift_ppe=noshift_ppe
        self.options.noshift_ppe2=noshift_ppe2
                
    def pre_recon_processing(self):
        if (not self.options.phasecorr_data):
            print('No Phasecorr Data') 
            return None   #something else#
        else:
            print('Processing FSE phase corrections...')
        #fetch param_dict from phase correction fid
        if (self.inputAcq.platform=="Varian"):
            phasecorrAcq = vrf.VarianAcquisition(self.options.phasecorr_data)
        else:
            print("Not ready for Bruker acquisition...")
        #fetch fid data from phase correction
        phase_fid_data = read_correction_fid(phasecorrAcq)

        #specify phase correction petable
        if (not self.options.phasecorr_table):
            print('No phase correction petable specified')
            default_table = os.path.join(self.options.phasecorr_data,
                                         'phasecorrtable_nf%d_ni%d'%(vrf.get_dict_value(phasecorrAcq.param_dict,'nf',6),
                                                                     vrf.get_dict_value(phasecorrAcq.param_dict,'ni',504)))
            phasecorr_petable_name = vrf.get_dict_value(phasecorrAcq.param_dict,'petable',default_table)   #fetch default table
        else:
            phasecorr_petable_name = self.options.phasecorr_table
        #petable checks
        t1_array = vrf.parse_petable_file(phasecorr_petable_name,'t1')
        t2_array = vrf.parse_petable_file(phasecorr_petable_name,'t2')
        nfid = phasecorrAcq.nfid
        etl = int(vrf.get_dict_value(phasecorrAcq.param_dict,'etl',6))
        nf = phasecorrAcq.nf
        ni_perchan = phasecorrAcq.ni*phasecorrAcq.nfid
        if any((t1_array[1::etl]!=0)&(t2_array[1::etl]!=0)):
            print('PE table format appears to be wrong...should have zeros on echo 2!!')
        #parse zeros from table
        ifid,itrain = nonzero( all( (reshape(t1_array,(nfid,ni_perchan*(nf//etl)//nfid,nf))==0)*
                                    (reshape(t2_array,(nfid,ni_perchan*(nf//etl)//nfid,nf))==0) ,axis=-1) )
        #"moments" from table
        M = moment_terms_from_petable(t1_array,t2_array,etl)
        M.shape = (M.shape[0],nfid,ni_perchan*(nf//etl)//nfid,etl)
        phase_fid_data.shape=(nfid,phase_fid_data.shape[1],ni_perchan*(nf//etl)//nfid,etl,phase_fid_data.shape[-1])
        #generate read shifts
        for imouse in range(phase_fid_data.shape[1]):
            self.maxroind = int( argmax(abs(mean( phase_fid_data[ifid,imouse,itrain,1,:] ,axis=0))) )
        for imouse in range(phase_fid_data.shape[1]):
            if (not self.options.no_even_odd_echo_roshift):
                even_odd_pixel_shift = ROshift(mean(mean(phase_fid_data[ifid,imouse,itrain,0::2,:],axis=0),axis=0),
                                               mean(mean(phase_fid_data[ifid,imouse,itrain,1::2,:],axis=0),axis=0))
                phase_fid_data[:,imouse,:,0::2,:] = apply_RO_shift(phase_fid_data[:,imouse,:,0::2,:], even_odd_pixel_shift/2.0) 
                phase_fid_data[:,imouse,:,1::2,:] = apply_RO_shift(phase_fid_data[:,imouse,:,1::2,:],-even_odd_pixel_shift/2.0) 
                self.even_odd_pixel_shift[imouse] = even_odd_pixel_shift
            elif (imouse==0): print('Not performing even_odd_roshift')   
            if (not self.options.no_even_odd_echo_phaseshift):
                self.ephase_corr[imouse,:] = Echo_phasecorr(mean(phase_fid_data[ifid,imouse,itrain,:,:],axis=0),alpha=self.alpha)
            elif (imouse==0): print('Not performing even_odd_echo_phaseshift')    
            if (not self.options.no_moment_terms_from_petable):
                self.pMfit[imouse,:] = Phase_shift_cal(phase_fid_data[:,imouse,:,:,:],M,ifid,itrain)
            elif (imouse==0): print('Not calculating moment_terms_from_petable')

    def gen_kspace(self,imouse=0):
        #grabs some useful numbers
        etl = int(vrf.get_dict_value(self.inputAcq.param_dict,'etl',6))
        nf = self.inputAcq.nf
        ni_perchan = self.inputAcq.ni*self.inputAcq.nfid
        nro = int(vrf.get_dict_value(self.inputAcq.param_dict,'np',630*2)/2)
        nrcvrs = self.inputAcq.nrcvrs
        t1_array = vrf.parse_petable_file(self.petable_name,'t1')
        t2_array = vrf.parse_petable_file(self.petable_name,'t2')
        #generate image moments to multiply with pMfit values
        M2 = moment_terms_from_petable(t1_array,t2_array,etl)
        #read in image k-space data
        kspace = zeros((nf*ni_perchan//etl,etl,nro),complex)
        for k in range(nf*ni_perchan//etl):
            fid_data,data_error = self.inputAcq.getdatafids(k*etl,(k+1)*etl,rcvrnum=imouse)
            kspace[k,:,:] = fid_data.copy()
        #apply read correction
        roramp = exp(1.j*2*pi*0.5*self.even_odd_pixel_shift[imouse]*(append(arange(nro/2),arange(-nro/2,0,1)))/nro)
        for k in range(kspace.shape[1]):
            kspace[:,k,:] = ifft(roramp**((-1)**(k%2))*self.ephase_corr[imouse,k]*
                                    fft(kspace[:,k,:],axis=-1),axis=-1)
        #apply phase correctionmaxroind
        phasecorr = exp(-1.j*sum(self.pMfit[imouse,:,newaxis,newaxis]*M2,axis=0))
        kspace = phasecorr[:,:,newaxis]*kspace 
        kspace1,kspace2 = split_kspace(self.inputAcq,kspace,self.petable_name,imouse,self.options)
        del kspace                       
        if (not self.options.no_Echo_shift_apply):
            kspace1,kspace2 = Echo_shift_apply(kspace1,kspace2,self.petable_name,self.inputAcq,self.maxroind)
        else: print('Not performing Echo_shift')  
        #average two acquisitions
        self.kspace = 0.5*(kspace1 + kspace2)
        del kspace1
        del kspace2
        #filter if requested
        if self.options.fermi_ellipse:
            self.kspace = fermi_ellipse_filter(self.kspace)
        elif self.options.fermi_pecirc:
            self.kspace = fermi_pecircle_filter(self.kspace)

 
#__________________________________________FUNCTIONS______________________________________#


###################### phase correction from short phasecorr experiment ########################################

def moment_terms_from_petable(t1_array,t2_array,etl):   
    trs_per_table = len(t1_array)//etl
    M = zeros((5,trs_per_table,etl),int)
    for k in range(trs_per_table): 
        cph = array((0,0,0,0,0),int)
        for j in range(etl):
            cph = -1*cph
            M[0:2,k,j] = cph[0:2]+array((t2_array[k*etl+j],t1_array[k*etl+j]),int)
            M[2:4,k,j] = cph[2:4]+array((t2_array[k*etl+j]**2,t1_array[k*etl+j]**2),int)
            M[4,k,j] = cph[4]+t2_array[k*etl+j]*t1_array[k*etl+j]
            cph[0:2] = cph[0:2]-array((t2_array[k*etl+j],t1_array[k*etl+j]),int)
            cph[2:4] = cph[2:4]+array((t2_array[k*etl+j]**2,t1_array[k*etl+j]**2),int)
            cph[4] = cph[4]+t2_array[k*etl+j]*t1_array[k*etl+j]
    return M    
    

#  get procpar parameters from phase correction petable        
def read_correction_fid(inputAcq):
    nro = int(vrf.get_dict_value(inputAcq.param_dict,'np',1)/2)
    etl = int(vrf.get_dict_value(inputAcq.param_dict,'etl',1))
    ni_perchan = inputAcq.ni*inputAcq.nfid
    raw_data = zeros((inputAcq.nfid,inputAcq.nmice,ni_perchan*(inputAcq.nf/etl)/inputAcq.nfid,inputAcq.nf,nro),complex) 
    for imouse in range (inputAcq.nmice):
        for j in range((ni_perchan*(inputAcq.nf//etl))):
            fid_data,data_error = inputAcq.getdatafids(j*inputAcq.nf,(j+1)*inputAcq.nf,rcvrnum=imouse)
            raw_data[j//(ni_perchan*(inputAcq.nf//etl)//inputAcq.nfid),imouse,j%(ni_perchan*(inputAcq.nf//etl)//inputAcq.nfid),:,:] = fid_data[:,:].copy()
    return raw_data 
    
   
def fitdiff(x,D_ph,dM):
    phfit = sum(x[:,newaxis]*dM,axis=0)
    return ravel(D_ph-phfit)    
    
#function to find RO shift between otherwise equivalent k-space lines
def ROshift(kline1,kline2,start=-2,stop=2,step=0.1):
    #phaseramp=exp(1.j*(arange(kline1.shape[-1])-kline1.shape[-1]/2)*2*pi/kline1.shape[-1])     
    phaseramp=exp(1.j*2*pi*(append(arange(kline1.shape[-1]/2),arange(-kline1.shape[-1]/2,0,1)))/kline1.shape[-1])
    steps=arange(start,stop,step)   
    C = zeros((len(steps),),float)
    for j in range(len(steps)):
        pixshift = steps[j]
        kline1adj = ifft((phaseramp)**(pixshift/2.0)*fft(kline1,axis=-1),axis=-1)    
        kline2adj = ifft((phaseramp)**(-pixshift/2.0)*fft(kline2,axis=-1),axis=-1)       
        C[j] = sum(abs(kline1adj)*abs(kline2adj))/(0.5*sum(abs(kline1adj)**2)+0.5*sum(abs(kline2adj)**2))
    i1 = argmax(C)
    istart = [0,i1-5][i1-5>0]
    iend = [len(steps),i1+5][i1+5<len(steps)]
    pfit=polyfit(steps[istart:iend],C[istart:iend],2)
    bestshift=-1*pfit[1]/(2.0*pfit[0])
    return bestshift

   
#readout shift calibration
def apply_RO_shift(kline,pixel_shift):               
    nro=kline.shape[-1]
    roramp = exp(1.j*2*pi*0.5*pixel_shift*(append(arange(nro//2),arange(-nro/2,0,1)))/nro)
    # second: roramp = exp(1.j*arange(nro)*2*pi*pixel_shift/nro) 
    # first: roramp = exp(1.j*(arange(nro)-nro/2)*2*pi*pixel_shift/nro) 
    klinemod = ifft(((roramp)*fft(kline,axis=-1)),axis=-1)
    #klinemod = fftshift(ifft(fftshift(roramp*fftshift(fft(fftshift(kline,axes=(-1,)),axis=-1),axes=(-1,)),axes=(-1,)),axis=-1),axes=(-1,))  
    return klinemod
   

    
def Echo_phasecorr(fid_echo_train,alpha=0.02):           
    maxroind = int( argmax(abs(fid_echo_train[1,:])) )
    ephase = unwrap(angle(fid_echo_train[:,maxroind]),axis=0)
    ephase_corr = exp(-1.j*ephase)
    eampl = abs(fid_echo_train[:,maxroind])
    eampl_corr = eampl/eampl[0]
    eampl_corr = eampl_corr/(eampl_corr**2+alpha)
    ephase_corr *= eampl_corr
    return ephase_corr

    
    
def Phase_shift_cal(imouse_phase_fid_data,M,ifid,itrain):         #imouse_phase_fid_data shape: (nfid,ni_perchan,nf,nro) 
    maxroind = int( argmax(abs(mean( imouse_phase_fid_data[ifid,itrain,1,:] ,axis=0))) ) 
    D_ph = median( angle((imouse_phase_fid_data[:,:,1::,maxroind-2:maxroind+3])*
                          conj(mean(imouse_phase_fid_data[ifid,itrain,1::,maxroind-2:maxroind+3],axis=0))),axis=-1)
    for j in range(D_ph.shape[0]):
        D_ph[j,itrain[j]::,:] = unwrap(D_ph[ifid[j],itrain[j]::,:],axis=-2)
        D_ph[j,0:itrain[j],:] = unwrap(D_ph[ifid[j],0:itrain[j],:][::-1,:],axis=-2)[::-1,:]
    pfit,cov_x,infodict,mesg,ier = leastsq(fitdiff,zeros(M.shape[0]),
                                            args=(reshape(D_ph[:,:,:],(prod(D_ph.shape),)),
                                                  reshape(M[:,:,:,1::],(M.shape[0],prod(M[:,:,:,1::].shape[1::])))
                                             ),full_output=True)
    return pfit           

  
###define kspace_2 planes   
def split_kspace(inputAcq,kspace,petable_name,imouse,options): 
    nv = int(vrf.get_dict_value(inputAcq.param_dict,'nv',504))
    nv2 = int(vrf.get_dict_value(inputAcq.param_dict,'nv2',504))
    t1 = vrf.parse_petable_file(petable_name,'t1')
    t2 = vrf.parse_petable_file(petable_name,'t2')
    kspace1 = petable_orderedpair_reordering(kspace[0:kspace.shape[-3]//2,:,:],('t1','t2'),petable_name=petable_name,matrix=(nv,nv2),
                                             index_start=0,index_end=(kspace.shape[-2]*kspace.shape[-3]//2))
    kspace1 = fov_adjustment(kspace1,options,inputAcq,imouse)
    print("Not performing ppe FOV adjustment")
    kspace2 = petable_orderedpair_reordering(kspace[kspace.shape[-3]//2:kspace.shape[-3],:,:],('t1','t2'),petable_name=petable_name,matrix=(nv,nv2),
                                             index_start=(kspace.shape[-2]*(kspace.shape[-3]//2)),index_end=(kspace.shape[-2]*kspace.shape[-3]))
    kspace2 = fov_adjustment(kspace2,options,inputAcq,imouse)
    return kspace1,kspace2    

###function to perform echo max shift    
def Echo_shift(kplane1,kplane2,minstep=-2,maxstep=2,stepsize=0.05):
    phaseramp1=exp(-1.j*2*pi*(append(arange(kplane1.shape[-1]/2),arange(-kplane1.shape[-1]/2,0,1)))/kplane1.shape[-1])
    phaseramp2=exp(-1.j*2*pi*(append(arange(kplane1.shape[-2]/2),arange(-kplane1.shape[-2]/2,0,1)))/kplane1.shape[-2])
    steps=arange(minstep,maxstep,stepsize) 
    xyvals = zeros((len(steps),len(steps),2),float)
    xyvals[:,:,1]=steps[newaxis,:]
    xyvals[:,:,0]=steps[:,newaxis]
    D = zeros((len(steps),len(steps)),float)
    for j in range(len(steps)):
        pixshift1 = steps[j]
        for k in range(len(steps)):
            pixshift2 = steps[k]          
            kplane1adj = ifft2((phaseramp2[:,newaxis])**(pixshift2/2.0)*
                               (phaseramp1[newaxis,:])**(pixshift1/2.0)*
                               fft2d(kplane1))
            kplane2adj = ifft2((phaseramp2[:,newaxis])**(-pixshift2/2.0)*
                               (phaseramp1[newaxis,:])**(-pixshift1/2.0)*
                               fft2d(kplane2))
            D[k,j] = sum(abs(kplane1adj)*abs(kplane2adj))/(0.5*sum(abs(kplane1adj)**2)+0.5*sum(abs(kplane2adj)**2))
    i_0 = argmax(D.flat)
    i_2,i_1 = (i_0/D.shape[1],i_0%D.shape[1])
    delta=2
    i_2start=[0,i_2-delta][i_2>=delta]
    i_2end = [len(steps)-1,i_2+delta][i_2<len(steps)-delta]
    i_1start=[0,i_1-delta][i_1>=delta]
    i_1end = [len(steps)-1,i_1+delta][i_1<len(steps)-delta]
    A = zeros(((i_2end-i_2start)*(i_1end-i_1start),6),float)
    A[:,0]=ravel(xyvals[i_2start:i_2end,i_1start:i_1end,1]**2)
    A[:,1]=ravel(xyvals[i_2start:i_2end,i_1start:i_1end,1]*xyvals[i_2start:i_2end,i_1start:i_1end,0])
    A[:,2]=ravel(xyvals[i_2start:i_2end,i_1start:i_1end,0]**2)
    A[:,3]=ravel(xyvals[i_2start:i_2end,i_1start:i_1end,1])
    A[:,4]=ravel(xyvals[i_2start:i_2end,i_1start:i_1end,0])
    A[:,5]=ones(A.shape[0])
    N = ravel(D[i_2start:i_2end,i_1start:i_1end])
    cfit,resids,rank,s=lstsq(A,N)
    maxx=(((cfit[1])*(cfit[4]))-(2*(cfit[3])*(cfit[2]))/((4*cfit[0]*cfit[2])-(cfit[1])**2))   # x= first from when dD/dx=0
    maxy=(((cfit[1]*cfit[3])-(2*cfit[4]*cfit[0]))/((4*cfit[0]*cfit[2])-((cfit[1])**2)))       # y= first from when dD/dy=0
    bestshift_2=array([maxy,maxx],float)    
    return bestshift_2    
    
#echo shift for centre only    
def Echo_shift_apply(kspace1,kspace2,petable_name,inputAcq,maxroind):   
    t1_array_petable = vrf.parse_petable_file(petable_name,'t1')
    t2_array_petable = vrf.parse_petable_file(petable_name,'t2')
    nv = kspace1.shape[-2] 
    nv2 = kspace1.shape[-3] 
    etl = int(vrf.get_dict_value(inputAcq.param_dict,'etl',6))
    nfid = inputAcq.nfid
    ni_perchan = inputAcq.ni*inputAcq.nfid
    nf = inputAcq.nf
    #defining petable1 petable2
    index1 = 0
    index2 = t1_array_petable.shape[0]//2
    index3 = t1_array_petable.shape[0]
    t1_array_petable1 = t1_array_petable[index1:index2]
    t1_array_petable2 = t1_array_petable[index2:index3]
    t2_array_petable1 = t2_array_petable[index1:index2]
    t2_array_petable2 = t2_array_petable[index2:index3]
    petable1_pemat=zeros((nv2,nv),int)
    petable1_pemat[t2_array_petable1+nv2//2-1,t1_array_petable1+nv//2-1]=1+arange(len(t1_array_petable1))%etl
    petable2_pemat=zeros((nv2,nv),int)
    petable2_pemat[t2_array_petable2+nv2//2-1,t1_array_petable2+nv//2-1]=1+arange(len(t1_array_petable2))%etl
    del t1_array_petable1
    del t1_array_petable2
    del t2_array_petable1
    del t2_array_petable2
    ROI1_centre = (petable1_pemat==petable1_pemat[nv2//2-1,nv//2-1])
    ROI2_centre = (petable2_pemat==petable2_pemat[nv2//2-1,nv//2-1])
    kspace1_centre = ROI1_centre[:,:,newaxis]*kspace1
    kspace2_centre = ROI2_centre[:,:,newaxis]*kspace2
    del ROI1_centre
    del ROI2_centre
    kspace1 = kspace1-kspace1_centre
    kspace2 = kspace2-kspace2_centre
    eopixshift_2 = Echo_shift((mean(kspace1_centre[:,:,maxroind-5:maxroind+5],axis=-1)),    
                     (mean(kspace2_centre[:,:,maxroind-5:maxroind+5],axis=-1)))
    eoramp_2 = exp(-1.j*2*pi*0.5*eopixshift_2[0]*append(arange(nv2//2),arange(-nv2//2,0,1))//nv2)[:,newaxis]* \
               exp(-1.j*2*pi*0.5*eopixshift_2[1]*append(arange(nv//2),arange(-nv//2,0,1))//nv)[newaxis,:]
    for j in range(kspace1_centre.shape[2]):
        kspace1_centre[:,:,j] = ifft2(eoramp_2*fft2(kspace1_centre[:,:,j]))
        kspace2_centre[:,:,j] = ifft2(eoramp_2**(-1)*fft2(kspace2_centre[:,:,j])) 
    kspace1 = kspace1 + kspace1_centre
    kspace2 = kspace2 + kspace2_centre
    del kspace1_centre
    del kspace2_centre
    return kspace1,kspace2
    

###################################### image recon ##############################################################    
    
# Example command:

#/home/bjnieman/source/vnmr/varian_recon_4/varian_recon.py /micehome/leigh/PythonPlay/cyl_fse_recon/Brian/fse3dmice_recon \
#              --phasecorr_data /projects/souris/leigh/fid/03feb14.fid_75um/phasecalibration \
#              --phasecorr_table /projects/souris/leigh/fid/03feb14.fid_75um/phasecalibration/phasecorr_272_nf8_ni272 \
#              --petable /projects/souris/leigh/fid/03feb14.fid_75um/LSN_E3_E4_pseudocyl_etl8_272_nf64_ni80 \
#              --no_even_odd_echo_roshift \
#              --no_even_odd_echo_phaseshift \
#              --no_moment_terms_from_petable \
#              --no_Echo_shift_apply \
#              --mouse_list 6 \
#              --no_fft \
#              /projects/souris/leigh/fid/03feb14.fid_75um\
#              /projects/souris/leigh/fid/03feb14.fid_75um/img_NEWphasecorr.mnc   
  


   

    

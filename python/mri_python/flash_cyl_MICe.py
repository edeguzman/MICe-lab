
from optparse import OptionGroup
import mri_python.recon_genfunctions as rgf
from mri_python.varian_read_file import parse_petable_file
from numpy import * #nonzero,empty,append,zeros,exp,arange,abs
from scipy.optimize import leastsq

def seq_specific_options(parser):
    optgroup = OptionGroup(parser,"Default Sequence","sequence-specific reconstruction options")
    optgroup.add_option("--petable_ordered_pairs", action="store_true", dest="petable_ordered_pairs",
                        default=0, help="phase encodes specified by (t1,t2) ordered pairs in petable")
    optgroup.add_option("--apowidth",type="float",dest="apowidth",
                        default=0.8, help="apodization transition width for partial fourier recon")    
    optgroup.add_option("--outputreps", action="store_true", dest="outputreps",
                        default=0, help="output repetitions rather than averaging")
    optgroup.add_option("--complexavg", action="store_true", dest="complexavg",
                        default=0, help="complex average of repetitions rather than magnitude")
    #optgroup.add_option("--petable_pe1", action="store_true", dest="petable_pe1",
    #                   default=0, help="use PE1 encoding order specified by t1 array in petable")
    #optgroup.add_option("--petable_pe2", action="store_true", dest="petable_pe2",
    #                    default=0, help="use PE2 encoding order specified by t1 array in petable")
    parser.add_option_group(optgroup)


class seq_reconstruction():
    def __init__(self,inputAcq,options,outputfile):
        self.options=options
        self.inputAcq=inputAcq
        self.outputfile=outputfile
        if ((self.options.petable=='') and (inputAcq.platform=="Varian")):
            self.petable_name=rgf.get_dict_value(inputAcq.param_dict,'petable','petable')
        elif ((self.options.petable=='') and (inputAcq.platform=="Bruker")):
            self.petable_name=rgf.get_dict_value(inputAcq.method_param_dict,'MICe_InputPetablePath','petable')
        else:
            self.petable_name=self.options.petable
        self.kspace=None
        self.image_data=None
        self.Pftacq=False
    def pre_recon_processing(self):
        pass
    def gen_kspace(self,imouse=0):
        rcvrelems = nonzero(self.inputAcq.rcvrmouse_mapping==imouse)[0]
        if (len(rcvrelems)==0):
            print("Requested reconstruction for mouse %d, which is not in receiver map"%imouse, self.inputAcq.rcvrmouse_mapping)
            raise SystemExit
        nTRs = len(parse_petable_file(self.petable_name,'t1'))
        self.inputAcq.nf = nTRs
        self.inputAcq.ni_perchan = 1
        nreps = rgf.get_dict_value(self.inputAcq.method_param_dict,"PVM_NRepetitions",1)
        self.kspace = empty([nreps,len(rcvrelems),self.inputAcq.ni_perchan,self.inputAcq.nf,self.inputAcq.nro],complex)
        for j,ircvr in enumerate(rcvrelems):
            self.kspace[:,j,...]=gen_kspace_local(self.inputAcq,ircvr,nTRs=nTRs)
        #self.kspace.shape = tuple([x for x in self.kspace.shape if x>1])
        if self.options.petable_ordered_pairs:
            self.kspace = rgf.petable_orderedpair_reordering(self.kspace,petable_arrays=('t1','t2'),petable_name=self.petable_name,\
                                                          matrix=(self.inputAcq.npe,self.inputAcq.npe2))
            print(self.kspace.shape)
        elif self.options.petable_pe1:
            self.kspace = rgf.petable_reordering(self.kspace,axis=-2,petable_array='t1',petable_name=self.petable_name)
        elif self.options.petable_pe2:
            self.kspace = rgf.petable_reordering(self.kspace,axis=-3,petable_array='t2',petable_name=self.petable_name)
        self.kspace = rgf.fov_adjustment(self.kspace,self.options,self.inputAcq,imouse)
        if (rgf.get_dict_value(self.inputAcq.method_param_dict,'PVM_EncPft',[1.0,1.0,1.0])[0]>1.0005):
            #first: zero-pad k-space
            nrofull = rgf.get_dict_value(self.inputAcq.method_param_dict,'PVM_Matrix',[1,1,1])[0]
            nroacq = rgf.get_dict_value(self.inputAcq.method_param_dict,'PVM_EncMatrix',[1,1,1])[0]
            npad= nrofull-nroacq
            self.kspace = append(zeros(self.kspace.shape[0:-1]+(npad,),complex),self.kspace,axis=-1)
            self.inputAcq.data_shape[-1] = nrofull
            self.Pftacq = True
        if self.options.fermi_ellipse:
            self.kspace = rgf.fermi_ellipse_filter(self.kspace)
    def recon(self):
        if ((self.Pftacq==True) and (not self.options.nofft)):
            #recon centre of k-space (with apodization) and obtain normed conjugate
            print("Partial fourier recon in the read direction...")
            nrofull = rgf.get_dict_value(self.inputAcq.method_param_dict,'PVM_Matrix',[1,1,1])[0]
            nroacq = rgf.get_dict_value(self.inputAcq.method_param_dict,'PVM_EncMatrix',[1,1,1])[0]
            filt = 1.0/(1.0+exp(-(nroacq-arange(nrofull))/self.options.apowidth))
            k_sym = self.kspace * filt
            img_sym = rgf.recon_3d(k_sym)
            del k_sym
            Pmult = exp(-1.j*angle(img_sym)) #conjugate(img_sym)/where(abs(img_sym))
            #multiply kspace by window function
            Kweight = 1.0/(1.0+exp((nrofull-nroacq-arange(nrofull))/self.options.apowidth))+ \
                      1.0/(1.0+exp((nroacq-arange(nrofull))/self.options.apowidth))
            weighted_image=rgf.recon_3d(self.kspace*Kweight).astype(complex)
            self.image_data=(weighted_image*Pmult).real
        else:
            rgf.default_recon(self)
        if ((self.options.complexavg) and (len(self.inputAcq.data_shape)>4)):
            thresh = percentile(abs(self.image_data).flat,75.0)
            nstep = 8
            g_coord=mgrid[0:self.image_data.shape[-3]:nstep,0:self.image_data.shape[-2]:nstep,0:self.image_data.shape[-1]:nstep]
            pts_log = (abs(self.image_data[0,0,ravel(g_coord[0]),ravel(g_coord[1]),ravel(g_coord[2])])>thresh)* \
                      (abs(self.image_data[0,1,ravel(g_coord[0]),ravel(g_coord[1]),ravel(g_coord[2])])>thresh)
            Ainds=reshape(g_coord,(3,len(pts_log)))[:,pts_log]
            Ivals1=mean(self.image_data[:,0,Ainds[0,:],Ainds[1,:],Ainds[2,:]],axis=0)
            Ivals2=mean(self.image_data[:,1,Ainds[0,:],Ainds[1,:],Ainds[2,:]],axis=0)
            Cvals=Ivals1*conj(Ivals2)
            Elem_phase_diff = arctan2(Cvals.imag,Cvals.real)
            pfit=Phase_shift_cal(Elem_phase_diff,Ainds)
            fullphaseimg = fullphasemap(pfit,self.image_data.shape)
            self.image_data = self.image_data[:,0,:,:,:]+self.image_data[:,1,:,:,:]*exp(1.j*fullphaseimg)
            self.image_data.shape = (self.image_data.shape[0],1,self.image_data.shape[-3],self.image_data.shape[-2],self.image_data.shape[-1])
        if (not self.options.outputreps):
            print("Averaging repetitions...")
            self.image_data = mean(abs(self.image_data),axis=0)
                
#dedicated gen_kspace function proves simpler than using gen_kspace_simple
def gen_kspace_local(inputAcq,ircvr,nTRs=None):
    print("Reading k-space data...")
    nrcvrs = inputAcq.nrcvrs
    no_ro = inputAcq.nro
    if ((no_ro<=0) or (no_ro>inputAcq.data_shape[-1])): no_ro = inputAcq.data_shape[-1] #np
    #if (len(inputAcq.data_shape)>2) and (inputAcq.data_shape[-3]>0):
    #    npe2 = inputAcq.data_shape[-3]
    #else:
    #    npe2 = 1
    #npe1 = inputAcq.data_shape[-2]
    if (inputAcq.platform=="Varian"):
        print("Not setup for Varian acquisition...")
        raise SystemExit
    elif (inputAcq.platform=="Bruker"): #force Bruker to artificial Varian format
        nf = nTRs 
        nreps = rgf.get_dict_value(inputAcq.method_param_dict,"PVM_NRepetitions",1)
    else:
        print("Input format not recognized.")
    raw_imgdata = zeros((nreps,1,nf,no_ro),complex)
    for k in range(nreps):
        print("%d / %d" % (k,nreps))
        fid_start = k*nf
        fid_end = (k+1)*nf
        fid_data,data_error = inputAcq.getdatafids(fid_start,fid_end,rcvrnum=ircvr)
        raw_imgdata[k,0,:,:] = fid_data[:,-no_ro:]
        if (data_error):
            data_fraction = float(data_error+fid_start)/float(nf*nreps)
            print('Error at %f%% through data set' % (100.0*data_fraction))
            break
    return raw_imgdata

def fitdiff(x,D_ph,M):
    phfit = sum(x[:,newaxis]*M,axis=0)
    return ravel(D_ph-phfit)    
    
def Phase_shift_cal(Phasediffvec,Vpos):         #imouse_phase_fid_data shape: (nfid,ni_perchan,nf,nro) 
    M=ones([10,Vpos.shape[1]])
    for j in range(3):
        M[j+1,:]=Vpos[j,:]
    for j in range(6):
        cinds=[[0,0],[1,1],[2,2],[0,1],[0,2],[1,2]][j]
        M[j+4,:]=Vpos[cinds[0],:]*Vpos[cinds[1],:]
    pfit,cov_x,infodict,mesg,ier = leastsq(fitdiff,zeros(M.shape[0]),
                                            args=(Phasediffvec,M),full_output=True)
    return pfit           

def fullphasemap(pfit,imgshape):
    g_coord=mgrid[0:imgshape[-3],0:imgshape[-2],0:imgshape[-1]]
    phasefit = pfit[0]+g_coord[0]*pfit[1]+g_coord[1]*pfit[2]+g_coord[2]*pfit[3]+ \
               g_coord[0]**2*pfit[4]+g_coord[1]**2*pfit[5]+g_coord[2]**2*pfit[6] + \
               g_coord[0]*g_coord[1]*pfit[7]+g_coord[0]*g_coord[2]*pfit[8]+g_coord[1]*g_coord[2]*pfit[9]
    return phasefit

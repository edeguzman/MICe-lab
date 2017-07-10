import numpy as N
from os import stat
from os import path as ospath
#from struct import calcsize as sizeof
from mri_python.varian_read_file import *
from mri_python.varian_fid_corrections import *
from pylab import hist,median
from scipy.signal import medfilt
from mri_python.recon_genfunctions import get_dict_value
from scipy.stats import mode
from scipy.interpolate import LSQUnivariateSpline
import collections

def get_sg_data(inputAcq,imouse):
    print("Retrieving self-gate data...")
    sgpts = get_dict_value(inputAcq.param_dict,'np',1)/2-get_dict_value(inputAcq.param_dict,'nro',1)
    if (sgpts<1): return 0
    sg_fids,data_error = inputAcq.getdatafids(0,product(inputAcq.data_shape[0:-1])/inputAcq.nrcvrs,rcvrnum=imouse,nrcvrs=inputAcq.nrcvrs,startpt=0,endpt=sgpts)
    return sg_fids

def self_resp_gate(navdata,tr,resp_dur):
    print("Generating respiratory gates from self-gated FIDs...")
    compi=N.array((0.+1.j),N.complex)
    sizes = N.shape(navdata)
    ref_base = N.zeros(sizes,N.complex)
    diff_navdata = N.zeros(sizes,N.complex)
    ref_width = 512
    n_fids_perblock = 128
    num_blocks = int(N.ceil(sizes[0]/float(n_fids_perblock)))
    for j in range(num_blocks):
        fid_start = j*n_fids_perblock
        fid_end = (j+1)*n_fids_perblock
        if (fid_end>sizes[0]):
            fid_end = sizes[0]
        start_pt = {0:0, 1:fid_start-(ref_width-n_fids_perblock)/2}[fid_start-(ref_width-n_fids_perblock)/2>0]
        end_pt = {0:sizes[0], 1:fid_end+(ref_width-n_fids_perblock)/2+1}[fid_end+(ref_width-n_fids_perblock)/2<sizes[0]-1]
        if start_pt==0:
            end_pt = ref_width
        if end_pt==sizes[0]:
            start_pt = sizes[0]-ref_width
        n_pts = end_pt-start_pt    
        reference_fid = N.sort( navdata[start_pt:end_pt,:].real, axis=0)[n_pts/4:3*n_pts/4+1:n_pts/4,:] + \
                compi * N.sort( navdata[start_pt:end_pt,:].imag, axis=0)[n_pts/4:3*n_pts/4+1:n_pts/4,:]
        ref_base[fid_start:fid_end,:] = reference_fid[1,:]
        diff_navdata[fid_start:fid_end,:] = (navdata[fid_start:fid_end,:] - reference_fid[1,:])/(reference_fid[2,:]-reference_fid[0,:])
    resp_sig = N.sum(abs(diff_navdata[:,0:5])/4.0,axis=-1)
    testfig = N.ravel(resp_sig)
    ref_val = median(testfig)
    close_testfigs = N.compress(N.less(testfig,20.0),testfig)
    numbins=400
    [binfreq,binloc,patches]=hist(close_testfigs,numbins)
    max_freq = N.maximum.reduce(binfreq)
    testfig_mode = binloc[ N.argmax(binfreq) ]
    quiet_cutoff = binloc[ N.maximum.reduce(N.nonzero( N.greater(binfreq,max_freq/8) )[0]) ] 
    if (quiet_cutoff<2*testfig_mode): 
        quiet_cutoff=2*testfig_mode
    print(quiet_cutoff)
    resp_gate_sig = N.less(resp_sig,quiet_cutoff)
    resp_gate_sig_orig = resp_gate_sig.copy()
    medfiltwindow = 2*(int(resp_dur/tr))+1
    resp_gate_sig = medfilt(resp_gate_sig_orig,medfiltwindow)
    return resp_gate_sig,resp_sig

def gen_kspace_sg(inputAcq,resp_gate_sig,resp_sig,mousenum,gate_resp=True):
    print("Generating gated k-space data...")
    no_ro = inputAcq.data_shape[-1]
    if (get_dict_value(inputAcq.param_dict,'sgflag','n')=='y'):
        no_ro = get_dict_value(inputAcq.param_dict,'nro',0)
    nechoes = get_dict_value(inputAcq.param_dict,'etl',1)
    nfid_pts = inputAcq.data_shape[-1]
    no_inner_loop = 1 #int(get_dict_value(method_param_dict,'$PVM_NMovieFrames',1))
    no_middle_loop = 1 #int(get_dict_value(image_param_dict,'no_middle_loop','1'))
    no_outer_loop = int(get_dict_value(inputAcq.param_dict,'nreps',1))
    nexpt = no_inner_loop*no_middle_loop*no_outer_loop
    if ( len(inputAcq.data_shape)>2 ):
        n_slowpe = inputAcq.data_shape[-3]
    else:
        n_slowpe = 1
    n_fastpe = inputAcq.data_shape[-2]
    if (nechoes>1):
        raw_imgdata = N.zeros((nechoes,n_slowpe,n_fastpe,no_ro),N.complex)
    else:
        raw_imgdata = N.zeros((n_slowpe,n_fastpe,no_ro),N.complex)
    traces_per_file = n_slowpe*n_fastpe*nexpt
    traces_per_slowpe = n_fastpe*no_inner_loop*no_middle_loop
    traces_per_fastpe = no_inner_loop
    #timing for echo locations
    sgpts = inputAcq.data_shape[-1]-get_dict_value(inputAcq.param_dict,'nro',1)
    start_pt = sgpts
    end_pt = inputAcq.data_shape[-1]
    nacq_pts = end_pt-start_pt
    if (nechoes>1): #need to get this setup so that it finds or reads pts automatically
        #start_pt=[194,895]; end_pt=[578,1279]; nacq_pts=[384,384];
        start_pt=[61,405]; end_pt=[317,661]; nacq_pts=[256,256]; #[68,426] [324,682] [256,256]
    nrcvrs = inputAcq.nrcvrs
    #start_pt = int(0.5 + start_time*sw*1e-3*ro_oversample)
    #nacq_pts = int(0.45+no_ro*fracacqwindow/imgacqwindow)
    #end_pt = int(start_pt + ro_oversample*nacq_pts) #(echotimes+acqwindow/2.0)*sw*1e-3*ro_oversample
    for j in range(n_slowpe):
        print("%d / %d" % (j,n_slowpe))
        fid_data = zeros((traces_per_slowpe*no_outer_loop,inputAcq.data_shape[-1]),complex)
        resp_flag = zeros((traces_per_slowpe*no_outer_loop,),float)
        resp_curr_sig = zeros((traces_per_slowpe*no_outer_loop,),float)
        for k in range(no_outer_loop): #range(no_outer_loop)
            fid_start = j*traces_per_slowpe+k*traces_per_slowpe*n_slowpe
            fid_end = (j+1)*traces_per_slowpe+k*traces_per_slowpe*n_slowpe
            fid_data_temp,data_error = inputAcq.getdatafids(fid_start,fid_end,rcvrnum=mousenum)
            if (data_error):
                data_fraction = float(j)/float(n_slowpe)
                print('Error at %f%% through data set' % (100.0*data_fraction))
                break
            fid_data[k*traces_per_slowpe:(k+1)*traces_per_slowpe,:]=fid_data_temp
            if (gate_resp):
                resp_flag[k*traces_per_slowpe:(k+1)*traces_per_slowpe] = resp_gate_sig[fid_start:fid_end].astype(N.float)
                resp_curr_sig[k*traces_per_slowpe:(k+1)*traces_per_slowpe] = resp_sig[fid_start:fid_end].astype(N.float)
            else:
                resp_flag[k*traces_per_slowpe:(k+1)*traces_per_slowpe] = N.ones((traces_per_slowpe,),N.float)
        fid_shape = (no_outer_loop*no_middle_loop,n_fastpe,traces_per_fastpe,inputAcq.data_shape[-1])
        resp_flag = N.reshape(resp_flag,fid_shape[:-1])
        resp_curr_sig = N.reshape(resp_curr_sig,fid_shape[:-1])
        fid_data = N.reshape(fid_data,fid_shape)
        norm_val = N.sum(N.sum(resp_flag,axis=0),axis=1)
        zero_inds = N.nonzero(norm_val==0)[0]
        for index in zero_inds:
            print("|")
            best_index = argmin(ravel(resp_curr_sig[:,index,:]))
            resp_flag[ best_index/traces_per_fastpe ,index, best_index%traces_per_fastpe ] = 1        
        norm_val = N.sum(N.sum(resp_flag,axis=0),axis=1)
        if (nechoes<=1):
            raw_imgdata[j,:,-nacq_pts:] = ( N.sum(N.sum( fid_data[:,:,:,start_pt:end_pt]* \
                                                         resp_flag[:,:,N.newaxis]       \
                                                        ,axis=0),axis=1)/norm_val[:,N.newaxis]).astype(N.complex)
        else:
            for k in range(nechoes):
                raw_imgdata[k,j,:,-nacq_pts[k]:] = ( N.sum(N.sum( fid_data[:,:,:,start_pt[k]:end_pt[k]]* \
                                                                  resp_flag[:,:,N.newaxis]       \
                                                                 ,axis=0),axis=1)/norm_val[:,N.newaxis]).astype(N.complex)

    return raw_imgdata

    

def gen_kspace_sg_orderedpairtable(inputAcq,resp_gate_sig,resp_sig,mousenum,gate_resp=True,
                                   petable=None,petable_arrays=('t1','t2'),grappafov=1,kpts_offset=0,outputreps=False,
                                   phasecorr=None,dcpl_info=None):
    print("Generating gated k-space data...")
    no_ro = inputAcq.data_shape[-1]
    if (get_dict_value(inputAcq.param_dict,'sgflag','n')=='y'):
        no_ro = get_dict_value(inputAcq.param_dict,'nro',0)
    nechoes = get_dict_value(inputAcq.param_dict,'etl',1)
    nfid_pts = inputAcq.data_shape[-1]
    nv = int( get_dict_value(inputAcq.param_dict,'nv',1) )
    nv2 = int( get_dict_value(inputAcq.param_dict,'nv2',1) )
    if (petable):
        t1array = parse_petable_file(petable,petable_arrays[0])
        t2array = parse_petable_file(petable,petable_arrays[1])
    else:
        t1array = (N.arange(nv*nv2)%nv)-nv//2+1
        t2array = (N.arange(nv*nv2)//nv)-nv2//2+1
    if (outputreps):
        gate_resp=False  #only makes sense to output reps without retrospective gating
        noutreps = int( mode(array( collections.Counter(t1array+nv//2-1+nv*(t2array+nv2//2-1)).most_common() )[:,1])[0][0] )  #count from table rather than looking in procpar
    else:
        noutreps = 1
    if (nechoes*noutreps>1):
        raw_imgdata = N.zeros((nechoes*noutreps,nv2,nv,no_ro),N.complex)
    else:
        raw_imgdata = N.zeros((nv2,nv,no_ro),N.complex)
    #timing for echo locations needs to be improved
    if (dcpl_info==None):
        start_pt = [inputAcq.data_shape[-1]-no_ro]
        end_pt = [inputAcq.data_shape[-1]]
    else:
        start_pt = [dcpl_info.rok0index-no_ro//2] #[0]
        end_pt= []
        for this_start_pt in start_pt:
            end_pt+=[this_start_pt + no_ro]

    nacq_pts = [end_pt[0]-start_pt[0]]
    if (nechoes>1): #need to get this setup so that it finds or reads pts automatically
        #start_pt=[194,895]; end_pt=[578,1279]; nacq_pts=[384,384];
        print('Warning: echo timing needs to be setup properly...')
        start_pt=[61,405]; end_pt=[317,661]; nacq_pts=[256,256]; #[68,426] [324,682] [256,256]
    nrcvrs = inputAcq.nrcvrs
    #start_pt = int(0.5 + start_time*sw*1e-3*ro_oversample)
    #nacq_pts = int(0.45+no_ro*fracacqwindow/imgacqwindow)
    #end_pt = int(start_pt + ro_oversample*nacq_pts) #(echotimes+acqwindow/2.0)*sw*1e-3*ro_oversample
    for j in range(nv2): 
        print("%d / %d" % (j,nv2))
        for k in range(nv): 
            fidinds = N.nonzero( (t2array==((j-nv2//2+1)*grappafov+kpts_offset))&
                                 (t1array==((k-nv//2+1)*grappafov+kpts_offset)) )[0]
            #print "gen_kspace_sg_orderedpairtable: ",len(fidinds)
            if (len(fidinds)==0):
                continue   
            nreps = len(fidinds)
            fid_data = zeros((nreps,inputAcq.data_shape[-1]),complex) #no_ro
            for q in range(nreps):
                #fid_data[q,:],data_error = get_vnmr_datafids(vnmrfidfilelist,fidinds[q],fidinds[q]+1,
                #                                             header_info,startpt=start_pt[0],endpt=end_pt[0],mouse_num=mousenum,nrcvrs=nrcvrs)
                fid_data[q,:] = get_corrected_datafids(inputAcq,fidinds[q],fidinds[q]+1,mouse_num=mousenum,phasecorr=phasecorr,dcpl_info=dcpl_info,dcpl_ppe_index=k-nv/2+1)                
            if (gate_resp):
                resp_flag = resp_gate_sig[fidinds].astype(N.int)
                resp_curr_sig = resp_sig[fidinds].astype(N.float)
            else:
                resp_flag = N.ones((nreps,),bool)
            norm_val = N.sum(resp_flag)
            if (norm_val==0):
                print("|")
                best_index = argmin(resp_curr_sig)
                resp_flag[ best_index ] = 1        
                norm_val = 1
            if (nechoes*noutreps<=1):
                raw_imgdata[j,k,-nacq_pts[0]:] = N.sum( fid_data[:,start_pt[0]:end_pt[0]]*resp_flag[:,N.newaxis]   \
                                                       ,axis=0)/norm_val
            else:
                repstep = nreps//noutreps
                rstarts = arange(noutreps)*repstep
                rends = (1+arange(noutreps))*repstep; rends=where(rends==rstarts,rstarts+1,rends); rends[-1]=len(fidinds)
                for cr in range(noutreps):
                    for q in range(nechoes):
                        raw_imgdata[cr*nechoes+q,j,k,-nacq_pts[q]:] = N.sum( fid_data[rstarts[cr]:rends[cr],
                                                                                      start_pt[q]:end_pt[q]]*resp_flag[rstarts[cr]:rends[cr],N.newaxis]       \
                                                                            ,axis=0)/N.sum(resp_flag[rstarts[cr]:rends[cr]])
    return raw_imgdata

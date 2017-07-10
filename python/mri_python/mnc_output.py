import numpy as N
#from py_minc import *
import imp
#from sys import path
#path.insert(0,"/axiom2/projects/software/arch/linux-precise/python/pyminc-0.41-py2.7.egg")
from pyminc.volumes.factory import *
from mri_python.recon_genfunctions import get_dict_value
from scipy.signal import medfilt
import string
import os

def Eulermat(psi,phi,theta):
    #definition of matrix here is odd, but consistent with the Agilent definition on our system, I think
    #need to verify for some orientations
    Emat = N.empty((3,3),float)
    Emat[0,0]=N.cos(psi)*N.sin(phi)-N.cos(theta)*N.cos(phi)*N.sin(psi)
    Emat[0,1]=-N.sin(psi)*N.sin(phi)-N.cos(theta)*N.cos(phi)*N.cos(psi)
    Emat[0,2]=N.cos(phi)*N.sin(theta)
    Emat[1,0]=-N.cos(psi)*N.cos(phi)-N.cos(theta)*N.sin(phi)*N.sin(psi)
    Emat[1,1]=N.sin(psi)*N.cos(phi)-N.cos(theta)*N.sin(phi)*N.cos(psi)
    Emat[1,2]=N.sin(phi)*N.sin(theta)
    Emat[2,0]=N.sin(theta)*N.sin(psi)
    Emat[2,1]=N.sin(theta)*N.cos(psi)
    Emat[2,2]=N.cos(theta)
    return Emat.T


    
#def add_scanparams_to_mncheader(image_param,filename):
#    print "Adding scan header to minc header..."
#    #full_scan_header = string.join([x.rstrip('\n').replace('"','') for x in image_param])
#    full_scan_header = string.join(image_param)
#    cmdstr = 'minc_modify_header -sinsert vnmr_scan_parameters:text=\"%s\" %s > /dev/null' \
#             % (full_scan_header,filename)
#    os.system(cmdstr)
   
def add_scanparams_to_mncheader(inputAcq,filename,imouse):
    print("Adding scan header to minc header...")
    nmice = inputAcq.nmice
    if (inputAcq.platform=="Varian"):
        dparamlist=['sw','tr','np','nv','nv2','nv_lores','nv2_lores','ni','nf','ne','etl','lro','lpe','lpe2',
                    'thk','ns','phi','psi','theta','fn','fn1','fn2','nmice','arrayelements','arraydim','nD','nleaf',
                    'nfid','npetables','frac_kx','esp','sfrq','te','ti','nt','ngainlo','gainlo','gainhi','gain',
                    'pss','mm_pro','mm_ppe','mm_ppe2','ppe','ppe2','pro']
        sparamlist=['array','seqfil','pslabel','orient','rfcoil','time_submitted','petable','rcvrs','seqcon','vrgain','acqtype',
                'username','datestr','mmconf','dp','console','console_rev','tn','protocols']
        image_param=inputAcq.param_dict
        platform_str="vnmr"
    elif (inputAcq.platform=="Bruker"):
        dparamlist=['PVM_EchoTime','PVM_RepetitionTime','PVM_NAverages','PVM_EffSWh','PVM_AcquisitionTime',
                    'PVM_SliceThick']
        sparamlist=['Method','OWNER','PVM_ScanTimeStr','PVM_Matrix','PVM_Fov','PVM_SPackArrSliceOrient','PVM_SPackArrReadOrient',
                    'PVM_SPackArrReadOffset','PVM_SPackArrPhase1Offset','PVM_SPackArrPhase2Offset','PVM_SPackArrSliceOffset',
                    'PVM_SPackArrSliceGap','PVM_SPackArrSliceDistance']
        image_param=inputAcq.method_param_dict
        platform_str="brkr"
    for ckey in dparamlist:
        cval = get_dict_value(image_param,ckey,0.0)
        if (hasattr(cval,'__len__')):
            if ((len(cval)==nmice) or (ckey=="mm_ppe") or (ckey=="mm_ppe2") or (ckey=="mm_pro")):
                cval = cval[imouse]
            else:
                cval = cval[0]
        cmdstr = 'minc_modify_header -dinsert %s:%s=%f %s > /dev/null' \
                 % (platform_str,ckey,cval,filename)
        os.system(cmdstr)
    cmdstr = 'minc_modify_header -dinsert %s:coil=%f %s > /dev/null' \
             % (platform_str,imouse+1,filename)
    os.system(cmdstr)
    for ckey in sparamlist:
        cval = str(get_dict_value(image_param,ckey,"NA"))
        cmdstr = 'minc_modify_header -sinsert %s:%s="%s" %s > /dev/null' \
                 % (platform_str,ckey,cval,filename)
        os.system(cmdstr)


def write_to_mnc_file(filename,data,inputAcq,options,manstarts=None,mansteps=None,imouse=0,mincversion=2):
    #output type
    nrcvrs_per_mouse=len(N.nonzero(inputAcq.rcvrmouse_mapping==imouse)[0])
    if options.real:
        if (data.dtype==complex):
            out_im=data.real
        else:
            out_im=data
        vType=[options.vType,"short"][options.vType==None]
    elif options.imag:
        out_im=data.imag; vType=[options.vType,"short"][options.vType==None]
    elif options.phase:
        out_im=N.arctan2(data.imag,data.real); vType=[options.vType,"short"][options.vType==None]
    elif (nrcvrs_per_mouse>1):
        print("Sum of squares reconstruction on multiple coils...")
        if ((inputAcq.nD==3) or ((inputAcq.nD==2) and (inputAcq.nslices>1))):
            rcvraxis=-4
        elif ((inputAcq.nD==2) and (inputAcq.nslices==1)):
            rcvraxis=-3
        elif (len(data.shape)>3):
            rcvraxis=-4
        else:
            rcvraxis=0
        if (len(data.shape)>inputAcq.nD):
            out_im=N.sqrt(N.sum(abs(data)**2,axis=rcvraxis))
        else: #assume only one coil has been received and that nD==ndims
            out_im=abs(data)
        vType=[options.vType,"ushort"][options.vType==None]
    else:
        out_im=abs(data); vType=[options.vType,"ushort"][options.vType==None]
    #dimensions (number and sizes)
    if (len(out_im.shape)>4):
        out_im.shape=(N.prod(out_im.shape[0:(len(out_im.shape)-3)]),out_im.shape[-3],out_im.shape[-2],out_im.shape[-1])
    dim_sizes=out_im.shape
    nD = len(dim_sizes)
    if (nD==2):
        nout_files=1
        dim_sizes = (1,dim_sizes[0],dim_sizes[1])
        out_im.shape = dim_sizes
    elif (nD<4):
        nout_files=1
    else:
        nout_files=dim_sizes[0]
        nD = 3
    if (inputAcq.platform=="Varian"):
        #fov and step size
        fov_array = [10.0*float(get_dict_value(inputAcq.param_dict,x,1.0)) for x in ['lpe2','lpe','lro']]
        (fov_ro,fov_pe,fov_pe2) = (fov_array[-1],fov_array[-2],fov_array[-3])
        #step & direction
        step_list=[fov_pe2/dim_sizes[-3],
                   fov_pe/dim_sizes[-2],
                   fov_ro/dim_sizes[-1],
                   1,1]
        [psi,phi,theta] = [get_dict_value(inputAcq.param_dict,x,0) for x in ['psi','phi','theta']]
        Emat = Eulermat(psi*N.pi/180,phi*N.pi/180,theta*N.pi/180) #did this end up Varian specific??
        dim_names=['zspace','yspace','xspace']    #default
        for j in range(3):
            dim_names[-1-j] = ['xspace','yspace','zspace'][N.argmax(abs(Emat[:,j]))]
            if (Emat[N.argmax(abs(Emat[:,j])),j]<0):
                step_list[2-j]=step_list[2-j]*(-1)
        #start point
        mm_pro = 10.0*get_dict_value(inputAcq.param_dict,'mm_pro',[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        mm_ppe = 10.0*get_dict_value(inputAcq.param_dict,'mm_ppe',[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        mm_ppe2 = 10.0*get_dict_value(inputAcq.param_dict,'mm_ppe2',[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        start_list=[N.sign(step_list[0])*mm_ppe2[imouse]-step_list[0]*(0.5*dim_sizes[-3]-0.5),
                    N.sign(step_list[1])*mm_ppe[imouse]-step_list[1]*(0.5*dim_sizes[-2]-0.5),
                    N.sign(step_list[2])*mm_pro[imouse]-step_list[2]*(0.5*dim_sizes[-1]-0.5),
                    0,0]        
    elif (inputAcq.platform=="Bruker"):
        #fov and step size
        fov_array = (get_dict_value(inputAcq.method_param_dict,"PVM_Fov",N.array([20.0,20.0,20.0],float)).astype(float))[::-1]
        (fov_ro,fov_pe) = (fov_array[-1],fov_array[-2])
        if (inputAcq.nD==3):
            fov_pe2 = fov_array[-3]
        elif (inputAcq.nD==2):
            fov_pe2 = get_dict_value(inputAcq.method_param_dict,"PVM_SPackArrSliceDistance",N.array([1.0,1.0,1.0],float))[0]*inputAcq.nslices
        #step & direction
        dim_names=['zspace','yspace','xspace']
        stepdirec=N.ones(3)
        Gmatrix_parse_row = lambda x: (['xspace','yspace','zspace'][N.argmax(abs(x))],[1,-1][x[N.argmax(abs(x))]<0])
        dim_names[0],stepdirec[0]=Gmatrix_parse_row(inputAcq.gmatrix[2])
        dim_names[1],stepdirec[1]=Gmatrix_parse_row(inputAcq.gmatrix[1])
        dim_names[2],stepdirec[2]=Gmatrix_parse_row(inputAcq.gmatrix[0])
        #dim_names=['zspace','yspace','xspace']    #default
        #sliceorient=get_dict_value(inputAcq.method_param_dict,'PVM_SPackArrSliceOrient',['axial'])[0]
        #dim_names[0]={'sagittal': 'xspace', 'coronal': 'yspace', 'axial': 'zspace'}[sliceorient]
        #readorient=get_dict_value(inputAcq.method_param_dict,'PVM_SPackArrReadOrient',['H_F'])[0]
        #dim_names[2]={'H_F': 'zspace','L_R': 'xspace', 'D_V': 'yspace'}[readorient] #what is yspace flag??? can we have F_H???
        #dim_names[1]=list(set(['zspace','yspace','xspace'])-set([dim_names[0],dim_names[2]]))[0]
        step_list=[fov_pe2*stepdirec[-3]/dim_sizes[-3],
                   fov_pe*stepdirec[-2]/dim_sizes[-2],
                   fov_ro*stepdirec[-1]/dim_sizes[-1],
                   1,1]
        #start point
        mm_pro = get_dict_value(inputAcq.method_param_dict,'MICe_ReadOffset',N.array([0.0,0.0,0.0,0.0],float))
        mm_ppe = get_dict_value(inputAcq.method_param_dict,'MICe_Phase1Offset',N.array([0.0,0.0,0.0,0.0],float))
        if (inputAcq.nD==3):
            #mm_ppe2 = get_dict_value(inputAcq.method_param_dict,'MICe_Phase2Offset',N.array([0.0,0.0,0.0,0.0],float))
            mm_ppe2 = get_dict_value(inputAcq.method_param_dict,'MICe_SliceOffset',N.array([0.0,0.0,0.0,0.0],float))
        elif (inputAcq.nD==2):
            mm_ppe2 = get_dict_value(inputAcq.method_param_dict,'MICe_SliceOffset',N.array([0.0,0.0,0.0,0.0],float))
        else:
            mm_ppe2 = N.array([0.0,0.0,0.0,0.0],float)
        start_list=[N.sign(step_list[0])*mm_ppe2[imouse]-step_list[0]*(0.5*dim_sizes[-3]-0.5),
                    N.sign(step_list[1])*mm_ppe[imouse]-step_list[1]*(0.5*dim_sizes[-2]-0.5),
                    N.sign(step_list[2])*mm_pro[imouse]-step_list[2]*(0.5*dim_sizes[-1]-0.5),
                    0,0]        
    else:
        print("Platform not recognized as Varian or Bruker")
    if (options.fft2d):
        #loop order: slice,pe1,ro
        if (inputAcq.platform=="Varian"):
            slice_sep = float(get_dict_value(inputAcq.param_dict,'thk',1.0)) #this isn't strictly the right parameter to use...
        elif (inputAcq.platform=="Bruker"):
            slice_sep = float(get_dict_value(inputAcq.method_param_dict,'PVM_SPackArrSliceDistance',[0.0])[0])
        start_list[0] = -dim_sizes[-3]*slice_sep/2.0 #this is not set up properly for 2D imaging yet...
        step_list[0] = slice_sep
    elif (options.nofft): #nofft
        if (fov_pe2>1e-6):
            steppe2=1.0/fov_pe2
            startpe2=-dim_sizes[-3]*steppe2/2
        else:
            steppe2=1.0
            startpe2=0
        step_list=[steppe2,1.0/fov_pe,1.0/fov_ro,1,1]
        start_list=[startpe2,-dim_sizes[-2]/fov_pe,-dim_sizes[-1]/fov_ro]
    #intensity scaling
    if (options.image_range_max-options.image_range_min>1e-4):
        print("Using user-input intensity range...")
        (Imin,Imax) = (options.image_range_min,options.image_range_max)
        for j in range(dim_sizes[0]):
            out_im[j] = N.where(N.greater(out_im[j],options.image_range_max),options.image_range_max,out_im[j]).astype(float)
            out_im[j] = N.where(N.less(out_im[j],options.image_range_min),options.image_range_min,out_im[j]).astype(float)
    else:
        N_voxels = N.size(out_im)
        flat_out = N.reshape(out_im,(N_voxels,))
        Imin=min(flat_out)
        if options.max_range:
            Imax=max(flat_out)
        else:
            plane_maxima = N.ravel(N.maximum.reduce(out_im,axis=-1))
            plane_maxima = medfilt(plane_maxima,11)
            Imax=1.1*N.maximum.reduce(plane_maxima)
            indices=N.nonzero(N.greater(flat_out,Imax))
            N.put(flat_out,indices,Imax)
    #override of starts and steps
    if (manstarts!=None):    
        start_list=list(manstarts)
    if (mansteps!=None):
        step_list=list(mansteps)
    #loop over output files
    for j in range(nout_files):
        #index output (for multiple echoes, repetitions, phases, etc)
        if (nout_files==1):
            outfile_name = filename
            imgshape = dim_sizes
        else:
            outfile_name = filename[:-4] + '_' + str(j) + '.mnc'
            imgshape = dim_sizes[1:]
        print("Writing to output file %s..." % outfile_name)
        #minc1 or minc2 output    
        if (mincversion!=2):
            print("Sorry, minc1 isn't an option.")

        spatDim = 3
        #print outfile_name,dim_names,imgshape,tuple(start_list[0:spatDim]),tuple(step_list[0:spatDim]),vType,out_im.dtype,out_im.shape
        mnc_vol = volumeFromDescription(outfile_name,dim_names,imgshape,tuple(start_list[0:spatDim]),tuple(step_list[0:spatDim]),volumeType=vType,dtype=out_im.dtype)
        if (nout_files==1):
            mnc_vol.setdata(out_im)
        else:
            mnc_vol.setdata(out_im[j])
        mnc_vol.writeFile()
        mnc_vol.closeVolume()
        #add scan params to mnc header
        add_scanparams_to_mncheader(inputAcq,outfile_name,imouse)
    return nout_files


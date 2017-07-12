# functions for opening, reading, closing Bruker fid files
import array as A
import os
from numpy import *
from struct import calcsize as sizeof
import re

program_name = 'bruker_read_file.py'

class FatalError(Exception):
    def __init__(self,args=None):
        self.msg = args

def get_dict_value(image_param_dict,key,default):
    retvalue=default
    if image_param_dict.has_key(key):
        retvalue=image_param_dict[key]
    return retvalue

class BrukerAcquisition():
    def __init__(self,inputpath,method_name="method",acqp_name="acqp",fid_name='fid'):
        self.platform="Bruker"
        #resolve file names and check existence
        if not (inputpath[-1]=='/'): inputpath = inputpath + '/'
        self.inputpath=inputpath
        method_file = os.path.join(inputpath,method_name)
        acq_file = os.path.join(inputpath,acqp_name)
        fid_file = os.path.join(inputpath,fid_name)
        try:
            if not ( os.path.exists(acq_file) ):
                raise FatalError("acqp file not found...")
            if not ( os.path.exists(fid_file) ):
                raise FatalError("fid file not found...")
            if not ( os.path.exists(method_file) ):
                raise FatalError("method file not found...")
        except FatalError as e:
            print('Error(%s):' % 'open_bruker_file', e.msg)
            raise SystemExit
        #generate bruker parameter dictionaries
        self.acq_param_dict=bruker_paramfile_to_dict(acq_file)
        self.method_param_dict=bruker_paramfile_to_dict(method_file)
        #specify datatype for "array" based read
        if (self.acq_param_dict['GO_raw_data_format']=='GO_16BIT_SGN_INT'): self.datatype = 'h'
        elif (self.acq_param_dict['GO_raw_data_format']=='GO_32BIT_FLOAT'): self.datatype = 'f'
        else: self.datatype = 'i' #default = 'GO_32BIT_SGN_INT'
        #determine data shape
        data_shape = ones(5,int) #(nreps*nframes*...),nrcvrs,nphase2|nslice,nphase,nro
        if ('PVM_EncMatrix' in self.method_param_dict):
            data_shape[-len(self.method_param_dict['PVM_EncMatrix'])::] = self.method_param_dict['PVM_EncMatrix'][::-1] #nphase2,nphase,nro
        elif ('PVM_SpecMatrix' in self.method_param_dict):
            #need to reshape this better
            data_shape[-3::] = array([1,1,self.method_param_dict['PVM_SpecMatrix'][0]],int) #leave 1's in for consistency
        if self.method_param_dict.has_key('PVM_EncActReceivers'):
           data_shape[-4] = self.method_param_dict["PVM_EncActReceivers"].count('On') #nrcvrs
        self.nrcvrs = data_shape[-4]
        nslices = sum(get_dict_value(self.method_param_dict,'PVM_SPackArrNSlices',[1]))
        if (nslices>1):
            data_shape[-3] *= nslices
        if self.method_param_dict.has_key('PVM_NRepetitions'):
            nreps = self.method_param_dict['PVM_NRepetitions']
            if (nreps>0): data_shape[-5] = nreps
        nframes = 0
        if self.method_param_dict.has_key('PVM_NMovieFrames'):
            if self.method_param_dict.has_key('PVM_MovieOnOff'):
                if (self.method_param_dict['PVM_MovieOnOff']=='On'):
                    nframes = self.method_param_dict['PVM_NMovieFrames']
        if (nframes>1):
            data_shape[-5] *= nframes
        #Do we need to verify the file size?
        #file_size = os.stat(fid_file)[6]
        #nfids = product(data_shape[:-1])
        #nropts = file_size/(2*nfids*sizeof(datatype))
        #if (data_shape[-1]!=nropts):
        #    print "Warning: inconsistency in file size and RO dimension (deriving number of RO points from file size)"
        #    data_shape[-1] = nropts
        self.data_shape=take(data_shape,nonzero(data_shape>1)[0])
        #for convenience, put some elements of data_shape into named variables
        self.nmice = get_dict_value(self.method_param_dict,'MICe_nmice',4)
        self.nro = data_shape[-1]
        self.npe = data_shape[-2]
        if (len(get_dict_value(self.method_param_dict,'PVM_EncMatrix',[1,1]))>2):
            self.npe2 = self.method_param_dict['PVM_EncMatrix'][2]
        else:
            self.npe2 = 1
        self.nslices = nslices
        #nD
        SpatDimEnum=get_dict_value(self.method_param_dict,'PVM_SpatDimEnum','<3D>')
        SpatDimEnumDict={'<3D>':3, '<2D>':2, '<1D>':1}
        self.nD=get_dict_value(SpatDimEnumDict,SpatDimEnum,3)
        #rcvrmouse_mapping: index is rcvr, value is mouse
        #need to read this from MICe_rcvrtomouse_mapping or similar
        complete_mapping=get_dict_value(self.method_param_dict,"MICe_receiver_to_coils",array([0,0,1,1,2,2,3,3],int))
        self.rcvrmouse_mapping=complete_mapping[ [i for i,onoff in enumerate(self.method_param_dict["PVM_EncActReceivers"]) if onoff=="On"] ] 
        #dummies for convenience when dealing with Varian
        self.nf=self.npe
        self.ni=[self.nslices,self.npe2][self.nD==3]
        self.ni_perchan=self.ni
        self.nfid=1 #may need to adjust this at some point
        #blocksize may be padded with zeroes depending on setting of GO_block_size
        base_blocksize=self.nrcvrs*self.nro*2*sizeof(self.datatype)
        if ((self.acq_param_dict["GO_block_size"]=="Standard_KBlock_Format") and not (fid_name[0:11]=="rawdata.job")): #temporary kluge, rawdata.jobN files don't follow KBlock format as fid files do
            self.blocksize=1024*( base_blocksize//1024 )
            if ((base_blocksize%1024)>0): self.blocksize+=1024
        else:
            self.blocksize=base_blocksize
        #define gradient matrix for orientation
        acqgradorient = get_dict_value(self.acq_param_dict,"ACQ_grad_matrix",None)
        acqpatientpos = get_dict_value(self.acq_param_dict,"ACQ_patient_pos",None)
        self.gmatrix = reshape(acqgradorient[0:9],(3,3))  #need to handle multiple orientations properly here
        if (acqpatientpos=="Head_Supine"):
            self.gmatrix[:,0]*=-1
            self.gmatrix[:,2]*=-1
        elif (acqpatientpos=="Head_Prone"):
            self.gmatrix[:,1]*=-1
            self.gmatrix[:,2]*=-1
        elif (acqpatientpos=="Head_Left"):
            self.gmatrix[:,2]*=-1
            self.gmatrix[:,[0,1]]=self.gmatrix[:,[1,0]] 
        elif (acqpatientpos=="Head_Right"):
            self.gmatrix[:,:]*=-1
            self.gmatrix[:,[0,1]]=self.gmatrix[:,[1,0]] 
        elif (acqpatientpos=="Foot_Supine"):
            self.gmatrix[:,:]=self.gmatrix[:,:] #no change for Foot_Supine
        elif (acqpatientpos=="Foot_Prone"):
            self.gmatrix[:,0]*=-1
            self.gmatrix[:,1]*=-1
        elif (acqpatientpos=="Foot_Left"):
            self.gmatrix[:,0]*=-1
            self.gmatrix[:,1]*=-1
            self.gmatrix[:,[0,1]]=self.gmatrix[:,[1,0]] 
        elif (acqpatientpos=="Foot_Right"):
            self.gmatrix[:,0]*=-1
            self.gmatrix[:,[0,1]]=self.gmatrix[:,[1,0]] 
        #Bruker convention seems to be to run -ve --> +ve on 1st phase encode,
        #  but then run +ve --> -ve on second phase encode; confusingly, this is
        #  sometimes recorded in the "ACQ_gradient_amplitude" value and sometimes in 
        #  the ppg instead (e.g., compare FLASH and RARE)  :(
        #  This may necessitate a sequence specific recon to set gmatrix properly
        #  on a case-by-case basis, matching the ACQ_gradient_amplitude, ppg, and petable
        #  statements, which may be switched in user sequences.
        #  Here, I am adding a -ve sign to the 2nd phase encode, so that the default 
        #  should handle the Bruker stock sequences properly
        self.gmatrix[2,:]*=-1
        print(self.gmatrix)
        #temporary handling of offsets
        default_hive_table=array([[-41,-41,0],[-41,41,0],[41,-40,0],[41,41,0]],float)
        if 1: #(not self.method_param_dict.has_key("MICe_ReadOffset")):
            self.method_param_dict["MICe_ReadOffset"]=array([sum(self.gmatrix[0,:]*default_hive_table[j]) for j in range(4)],float)
        if 1: #(not self.method_param_dict.has_key("MICe_Phase1Offset")):
            self.method_param_dict["MICe_Phase1Offset"]=array([sum(self.gmatrix[1,:]*default_hive_table[j]) for j in range(4)],float)
        if 1: #(not self.method_param_dict.has_key("MICe_Phase2Offset")):
            self.method_param_dict["MICe_Phase2Offset"]=array([sum(self.gmatrix[2,:]*default_hive_table[j]) for j in range(4)],float)
        if 1: #(not self.method_param_dict.has_key("MICe_SliceOffset")):
            self.method_param_dict["MICe_SliceOffset"]=array([sum(self.gmatrix[2,:]*default_hive_table[j]) for j in range(4)],float)
        #finally, open fid file for reading data
        self.fidfilehandle=open(fid_file,'r')
    def close(self):
        self.fidfilehandle.close()
    def getdatafids(self,fid_start,fid_end,rcvrnum=None):
        datasize = sizeof(self.datatype)
        nrcvrs = self.nrcvrs
        if (rcvrnum is None):
            crcvrs = arange(self.nrcvrs)
        else:
            try:
                test = len(rcvrnum)
                crcvrs = array(rcvrnum,int)
            except TypeError:
                crcvrs = array([rcvrnum],int)
        #complex data array
        complex_data = zeros((len(crcvrs),fid_end-fid_start,self.nro),complex)
        #read in data
        data_error = 0
        for k in range(len(crcvrs)):
            for j in range(fid_end-fid_start):
                self.fidfilehandle.seek(self.blocksize*(fid_start+j)+2*self.nro*datasize*crcvrs[k],0)
                bindata=A.array(self.datatype)
                try:
                    bindata.read(self.fidfilehandle,2*self.nro)
                except EOFError:
                    print('Error(%s): Missing data in file!' % program_name)
                    data_error = j
                    break
                complex_data[k,j,:]=array(bindata[0:2*self.nro:2],float)+1.j*array(bindata[1:2*self.nro+1:2],float)
        if (len(crcvrs)==1):
            complex_data.shape=(fid_end-fid_start,self.nro)
        return complex_data,data_error
        
def bruker_paramfile_to_dict(brukerparamfile):
    bpfh = open(brukerparamfile,'r')
    text_lines = bpfh.readlines()
    bpfh.close()  
    param_dict={}
    ntext = 0
    curr_line = 0
    nlines = len(text_lines)
    while (curr_line<nlines):
        line = text_lines[curr_line]
        if ((line[0:7]=='$$ @vis') or (line[0:7]=='$$ ')):
            curr_line += 1
            continue
        if line[0:3]=='##$':
            words = line[3:].split('=') 
            try:
                val = int(words[1])
            except ValueError as e:
                try:
                    val = float(words[1])
                except ValueError as f:
                    val = words[1][:-1]
                    if (len(val)==0):
                        curr_line += 1
                        continue
                    if (words[0]=="PVM_SliceGeo"): #handle SliceGeo as special case since it is important for orientation info
                        slicegeostr = ''                 #expand this to handle multiple slice orientations? currently just grabs first one...
                        while 1:
                            curr_line+=1
                            line = text_lines[curr_line]
                            if (line[0:2]=='##')|(line[0:2]=='$$'): break
                            slicegeostr += line[:-1] #exclude new line char at the end of the line
                        param_dict['slicegeostr']=slicegeostr
                        c_regexpr=re.compile("\(\(\((?P<Gmatrix>.*?),(?P<Offset>.*?)\),(?P<Fov>.*?), <(?P<Row1>.*?)> <(?P<Row2>.*?)> <(?P<Row3>.*?)>,(?P<ustr1>.*?)\),(?P<ustr2>.*?),(?P<ustr3>.*?),(?P<ustr4>.*?),(?P<rsltn>.*?),(?P<ustr5>.*?),(?P<ustr6>.*?)\)",re.DOTALL)
                        c_matchiter=c_regexpr.finditer(slicegeostr)
                        val=[]
                        for c_matches in c_matchiter:
                            valE={'Gmatrix':  reshape(array(map(float,c_matches.group('Gmatrix').split())),(3,3)),
                                 'Offset':   array(map(float,c_matches.group('Offset').split())),
                                 'Fov':      array(map(float,c_matches.group('Fov').split())),
                                 'RowOrder': [c_matches.group('Row1'),c_matches.group('Row2'),c_matches.group('Row3')],
                                 'ustr1':    c_matches.group('ustr1'),  #sort out these later if needed
                                 'ustr2':    c_matches.group('ustr2'),
                                 'ustr3':    c_matches.group('ustr3'),
                                 'ustr4':    c_matches.group('ustr4'),
                                 'rsltn':    c_matches.group('rsltn'),
                                 'ustr5':    c_matches.group('ustr5'),
                                 'ustr6':    c_matches.group('ustr6')}
                            GmatOrder = [[i for i,j in enumerate(valE["RowOrder"]) if 'read' in j][0],
                                         [i for i,j in enumerate(valE["RowOrder"]) if 'phase' in j][0],
                                         [i for i,j in enumerate(valE["RowOrder"]) if 'slice' in j][0]]
                            valE['GmatrixRPS'] = valE['Gmatrix'][GmatOrder,:]
                            val.append(valE)
                        curr_line -= 1
                    elif (val[0]=='(')&(val[-1]==')'):
                        elems = []
                        while 1:
                            curr_line+=1
                            line = text_lines[curr_line]
                            if (line[0:2]=='##')|(line[0:2]=='$$'): break
                            elem_words = line[:-1].split()
                            elems.extend(elem_words)
                        curr_line -= 1
                        try:
                            type_array = 'int'
                            for m in range(len(elems)):
                                val = int(elems[m])
                        except ValueError as e:
                            try:
                                type_array = 'float'
                                for m in range(len(elems)):
                                    val = float(elems[m])
                            except ValueError as f:
                                type_array = 'str'
                        if (type_array=='int'):
                            val = array([int(x) for x in elems],int)
                        elif (type_array=='float'):
                            val = array([float(x) for x in elems],float)
                        else:
                            val = elems
                    elif (val[0]=='('):
                        elems = val.split(',')
                        elems[0] = elems[0].strip('(')
                        while 1:
                            curr_line+=1
                            line = text_lines[curr_line]
                            if (line[0:2]=='##')|(line[0:2]=='$$'): break
                            elem_words = line[:-1].split()
                            elems.extend(elem_words)
                        elems[-1] = elems[-1].strip(')')
                        val = elems    
            param_dict[words[0]] = val
            curr_line += 1
            continue
        curr_line += 1
    return param_dict


#######################################################################################3
#keep the following two functions for compatibility with much older code
#shouldn't need to be used anymore though
#######################################################################################3

def open_bruker_file(inputpath):
    bacq = BrukerAcquisition(inputpath)
    return bacq.fidfilehandle,bacq.data_shape,bacq.datatype,bacq.acq_param_dict,"",bacq.method_param_dict,""

def get_bruker_datafids(brukerfile,fid_start,fid_end,nro,datatype):
    comp_i = array((0.+1.j),complex)
    fidarray = A.array(datatype)
    datasize = fidarray.itemsize
    #read in data to a data array
    brukerfile.seek(2*datasize*fid_start*nro,0)
    data_elems = nro*(fid_end-fid_start)
    complex_data = zeros((data_elems,),complex)
    data_error = 0
    for j in range(fid_end-fid_start):
        bindata=A.array(datatype)
        try:
            bindata.read(brukerfile,2*nro)
        except EOFError:
            print('Error(%s): Missing data in file!' % program_name)
            data_error = j
            break
        complex_data[j*nro:(j+1)*nro]=array(bindata[0:2*nro:2],float)+comp_i*array(bindata[1:2*nro+1:2],float)
    #force data to standard shape: reordering/sorting performed later
    data_shape = array((fid_end-fid_start,nro))
    complex_data = reshape(complex_data,tuple(data_shape))
    return complex_data,data_error

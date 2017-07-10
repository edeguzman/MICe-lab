
from optparse import OptionGroup
import mri_python.recon_genfunctions as rgf
from numpy import nonzero,empty

def seq_specific_options(parser):
    optgroup = OptionGroup(parser,"Default Sequence","sequence-specific reconstruction options")
    optgroup.add_option("--cen_out_pe1", action="store_true", dest="cen_out_pe1",
                        default=0, help="specify centre-out acq. order in first phase encode")
    optgroup.add_option("--cen_out_pe2", action="store_true", dest="cen_out_pe2",
                        default=0, help="specify centre-out acq. order in second phase encode")
    optgroup.add_option("--petable_pe1", action="store_true", dest="petable_pe1",
                       default=0, help="use PE1 encoding order specified by t1 array in petable")
    optgroup.add_option("--petable_pe2", action="store_true", dest="petable_pe2",
                        default=0, help="use PE2 encoding order specified by t1 array in petable")
    optgroup.add_option("--petable_ordered_pairs", action="store_true", dest="petable_ordered_pairs",
                        default=0, help="phase encodes specified by (t1,t2) ordered pairs in petable")
    optgroup.add_option("--dcshiftcorr",action="store_true",dest="dcshiftcorr",
                       default=0, help="attempt to remove DC artifact from images during recon")
    parser.add_option_group(optgroup)

    #parser.add_option("--avg_multi_fid",action="store_true",dest="avg_multi_fid",
    #                  default=0, help="average multiple fid files prior to reconstruction")


class seq_reconstruction():
    def __init__(self,inputAcq,options,outputfile):
        self.options=options
        self.inputAcq=inputAcq
        self.outputfile=outputfile
        if ((self.options.petable=='') and (inputAcq.platform=="Varian")):
            self.petable_name=rgf.get_dict_value(inputAcq.param_dict,'petable','petable')
        else:
            self.petable_name=self.options.petable
        self.kspace=None
        self.image_data=None
    def pre_recon_processing(self):
        pass
    def gen_kspace(self,imouse=0):
        rcvrelems = nonzero(self.inputAcq.rcvrmouse_mapping==imouse)[0]
        if (len(rcvrelems)==0):
            print("Requested reconstruction for mouse %d, which is not in receiver map"%imouse, self.inputAcq.rcvrmouse_mapping)
            raise SystemExit
        nreps = rgf.get_dict_value(self.inputAcq.method_param_dict,"PVM_NRepetitions",1)
        self.kspace = empty([len(rcvrelems),nreps*self.inputAcq.ni*self.inputAcq.nfid,self.inputAcq.nf,self.inputAcq.nro],complex)
        for j,ircvr in enumerate(rcvrelems):
            self.kspace[j] = rgf.gen_kspace_simple(self.inputAcq,ircvr)
        self.kspace.shape = tuple([x for x in self.kspace.shape if x>1])
        if self.options.cen_out_pe1:
            self.kspace=rgf.centre_out_reordering(self.kspace,-2)
        if self.options.cen_out_pe2:
            self.kspace=rgf.centre_out_reordering(self.kspace,-3)
        if self.options.petable_pe1:
            self.kspace = rgf.petable_reordering(self.kspace,axis=-2,petable_array='t1',petable_name=self.petable_name)
        if self.options.petable_pe2:
            self.kspace = rgf.petable_reordering(self.kspace,axis=-3,petable_array='t2',petable_name=self.petable_name)
        if self.options.petable_ordered_pairs:
            self.kspace = rgf.petable_orderedpair_reordering(self.kspace,petable_arrays=('t1','t2'),petable_name=self.petable_name,\
                                                          matrix=(self.inputAcq.npe,self.inputAcq.npe2))
        self.kspace = rgf.fov_adjustment(self.kspace,self.options,self.inputAcq,imouse)
        if self.options.dcshiftcorr: #do we use this?
            self.kspace,dcoff = rgf.DCartcorr(self.kspace,self.param_dict)
        if self.options.fermi_ellipse:
            self.kspace = rgf.fermi_ellipse_filter(self.kspace)
        elif self.options.fermi_pecirc:
            self.kspace = rgf.fermi_pecircle_filter(self.kspace) 



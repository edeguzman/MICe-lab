
from optparse import OptionGroup
import mri_python.recon_genfunctions as rgf
import mri_python.varian_fid_corrections as vfc
import mri_python.resp_sg as rsg
import mri_python.coil_decouple as cd

def seq_specific_options(parser):
    optgroup = OptionGroup(parser,"ge3dmice_sg_cylA","sequence-specific reconstruction options")
    optgroup.add_option("--petable_ordered_pairs", action="store_true", dest="petable_ordered_pairs",
                        default=0, help="phase encodes specified by (t1,t2) ordered pairs in petable")
    optgroup.add_option("--self_resp_gate",action="store_true",dest="self_resp_gate",
                       default=0, help="use FID nav for resp. gating")
    optgroup.add_option("--outputreps",action="store_true",dest="outputreps",
                       default=0, help="output all reps from a self-gated sequence separately") 
    optgroup.add_option("--phasedriftcorr",action="store_true",dest="phasedriftcorr",
                       default=0, help="correct phase drift based on repeated zero acqs in table") 
    optgroup.add_option("--rep_position_corr",action="store_true",dest="rep_position_corr",
                       default=0, help="correct phase encode shifts between multiple repetitions") 
    optgroup.add_option("--grappa_coil_decouple",action="store_true",dest="grappa_coil_decouple",
                       default=0, help="correct coupling between coils (i.e., in saddle coil array)") 
    optgroup.add_option("--grappa_coil_groupings",type="string",dest="grappa_coil_groupings",
                       default="0,2,5;1,3,4,6", help="groups of cross-talking coils (example: 0,2,5;1,3,4,6 )") 
    optgroup.add_option("--dcplppeadj",type="string",dest="dcplppeadj_string",
                       default=None, help="ppe position offset for decoupling mask")
    optgroup.add_option("--grappa_kpts_offset",type="int",dest="grappa_kpts_offset",
                       default=0, help="grappa to k-space grid offset (same in both directions: option to be deleted!)")
    optgroup.add_option("--dcpl_profs_output_name",type="string",dest="dcpl_profs_output_name",
                       default=None, help="filename for output of grappa profiles and mask contours (useful for debugging decoupling)")
    optgroup.add_option("--pre_applied_pe1_shift",action="store_true",dest="pre_applied_pe1_shift",
                       default=False, help="fid data has pre-applied pe1 shift")
    parser.add_option_group(optgroup)
    
    
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
        if (not self.options.dcplppeadj_string):
            self.dcplppeadj = None
        else:
            self.dcplppeadj = [float(x) for x in self.options.dcplppeadj_string.split(',')]

    def pre_recon_processing(self):
        if (self.options.phasedriftcorr):
            self.phasedriftcorr = vfc.phase_drift_corr(self.inputAcq,self.petable_name)
        else:
            self.phasedriftcorr = None
        if (self.options.rep_position_corr):
            self.phasedriftcorr = vfc.rep_pos_corr(self.inputAcq,self.petable_name,phasedriftcorr=self.phasedriftcorr)    
        if (self.options.grappa_coil_decouple):
            self.dcpl_info = cd.gen_decoupling_info(self.inputAcq,self.petable_name,
                                                    cplgrps_string=self.options.grappa_coil_groupings,
                                                    dcplppeadj=self.dcplppeadj,
                                                    dcpl_profs_output_name=self.options.dcpl_profs_output_name,
                                                    remove_ppeshift=self.options.pre_applied_pe1_shift)
        else:
            self.dcpl_info = None

    def gen_kspace(self,imouse=0):
        sgflag = rgf.get_dict_value(self.inputAcq.param_dict,'sgflag','n') 
        if (sgflag=='y') and self.options.self_resp_gate:
            sg_fids = rsg.get_sg_data(self.inputAcq,imouse)
            #from pylab import plot,show
            #plot(sg_fids[:,3].real)
            #plot(sg_fids[:,3].imag)
            #show()
            tr = rgf.get_dict_value(self.inputAcq.param_dict,'tr',0.05)
            resp_dur = 0.2
            resp_gate_sig,resp_sig = rsg.self_resp_gate(sg_fids,tr,resp_dur)
        else:
            resp_gate_sig=None
            resp_sig=None
        if (self.options.petable_ordered_pairs):
            grappafov=int(rgf.get_dict_value(self.inputAcq.param_dict,'grappafov',1))
            self.kspace = rsg.gen_kspace_sg_orderedpairtable(self.inputAcq,
                                                             resp_gate_sig,resp_sig,imouse,
                                                             gate_resp=self.options.self_resp_gate,
                                                             petable=self.petable_name,grappafov=grappafov,
                                                             kpts_offset=self.options.grappa_kpts_offset,
                                                             outputreps=self.options.outputreps,
                                                             phasecorr=self.phasedriftcorr,
                                                             dcpl_info=self.dcpl_info)
        elif (sgflag=='y'):
            self.kspace = rsg.gen_kspace_sg(self.inputAcq,resp_gate_sig,resp_sig,imouse)
        else: #sgflag='n' and no petable_ordered_pairs
            self.kspace = rgf.gen_kspace_simple(self.inputAcq,imouse)
        self.kspace = rgf.fov_adjustment(self.kspace,self.options,self.inputAcq,imouse)

        if self.options.fermi_ellipse:
            self.kspace = rgf.fermi_ellipse_filter(self.kspace)
        elif self.options.fermi_pecirc:
            self.kspace = rgf.fermi_pecircle_filter(self.kspace) 
    
    
    

        


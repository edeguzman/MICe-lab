#!/usr/bin/env python3
#
# mri_recon.py
#
# reconstruction of varian or bruker MR images with .mnc output
#
# Created June 2007


from optparse import OptionParser, Option, OptionValueError, OptionGroup
import importlib
import sys
import os
import re
import glob
import numpy as N
import zipimport
from mri_python.recon_genfunctions import FatalError,get_dict_value,default_recon
import mri_python.varian_read_file as vrf
import mri_python.bruker_read_file as brf
from mri_python.mnc_output import write_to_mnc_file

program_name = 'mri_recon.py'
DEFAULT_RECON_TYPE = 'seq_default'

############################################################################################################################################
    
def generate_option_parser_and_seq_module(recontypestring=None):
    usage = """%s <recon type flag> [options] <inputdirectory> <outputfile> 
   or  %s --help
   or  %s <recon type flag> --help

%s is a script for performing recon operations on varian fid files
(the input directory must contain both procpar and fid files)
"""
    usage = usage % ((program_name, )*4)
    parser = OptionParser(usage)
    #Generic recon options
    optgroup = OptionGroup(parser,"Generic","general info and reconstruction type options")
    optgroup.add_option("--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite existing output file")
    optgroup.add_option("--no_fft", action="store_true", dest="nofft",
                        default=0, help="output kspace data to file")
    optgroup.add_option("--1d_fft", action="store_true", dest="fft1d",
                        default=0, help="perform only 1d fft")
    optgroup.add_option("--2d_fft", action="store_true", dest="fft2d",
                        default=0, help="perform 2d fft reconstruction")
    optgroup.add_option("--3d_fft", action="store_true", dest="fft3d",
                       default=1, help="perform 3d fft reconstruction")
    optgroup.add_option("--large_data_recon", action="store_true", dest="large_data_recon",
                        default=0, help="perform 'in-place' 3dfft to save memory")
    optgroup.add_option("--real",action="store_true",dest="real",
                        default=0, help="write only real part of data to file")
    optgroup.add_option("--imag",action="store_true",dest="imag",
                        default=0, help="write only imag part of data to file")
    optgroup.add_option("--phase",action="store_true",dest="phase",
                        default=0, help="write phase of image to file")
    optgroup.add_option("--petable", type="string", dest="petable", metavar="petable",default='',
                        help="petable name and location if not specified by parameter files and current directory")
    optgroup.add_option("--mouse_list",type="str",dest="mouse_list",
                        default=None, help="reconstruct list of mice instead of all mice (zero-based indexing, comma-delimited, no spaces)")
    optgroup.add_option("--procpar_file_name", type="string", dest="procpar_file_name", metavar="procpar",default="procpar",
                        help="procpar file name for Varian acquisition (if not simply 'procpar')")
    optgroup.add_option("--method_file_name", type="string", dest="method_file_name", metavar="method",default="method",
                        help="method file name for Bruker acquisition (if not simply 'method')")
    optgroup.add_option("--acqp_file_name", type="string", dest="acqp_file_name", metavar="acqp",default="acqp",
                        help="acqp file name for Bruker acquisition (if not simply 'acqp')")
    optgroup.add_option("--fid_file_name", type="string", dest="fid_file_name", metavar="fid",default="fid",
                        help="fid file name for Bruker acquisition (if not simply 'fid')")                        
    parser.add_option_group(optgroup)
    #FOV options
    optgroup = OptionGroup(parser,"FOV","field-of-view shift options")    
    optgroup.add_option("--noshift_ppe",action="store_true",dest="noshift_ppe",
                        default=0, help="do not perform ppe shift described by mm_ppe (default: perform shift)")
    optgroup.add_option("--noshift_ppe2",action="store_true",dest="noshift_ppe2",
                        default=0, help="do not perform ppe shift described by mm_ppe (default: perform shift)")
    optgroup.add_option("--fov_shift_ro",type="string",dest="fov_shift_ro",
                        default="0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0", help="field-of-view shift in RO direction (# pixels, may be fractional, comma-separated list with no spaces)")
    optgroup.add_option("--fov_shift_pe1",type="string",dest="fov_shift_pe1",
                       default="0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0", help="field-of-view shift in PE1 direction (# pixels, may be fractional, comma-separated list with no spaces)")
    optgroup.add_option("--fov_shift_pe2",type="string",dest="fov_shift_pe2",
                      default="0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0", help="field-of-view shift in PE2 direction (# pixels, may be fractional, comma-separated list with no spaces)")
    parser.add_option_group(optgroup)
    #Filters                  
    optgroup = OptionGroup(parser,"Filters","filter options")    
    optgroup.add_option("--fermi_ellipse",action="store_true",dest="fermi_ellipse",
                       default=0, help="apply fermi ellipse filter")
    optgroup.add_option("--fermi_pecirc",action="store_true",dest="fermi_pecirc",
                       default=0, help="apply fermi circle filter to PE1 and PE2")
    parser.add_option_group(optgroup)
    #Common output options
    optgroup = OptionGroup(parser,"Output","common output options")
    optgroup.add_option("--max_range",action="store_true",dest="max_range",
                       default=0, help="set output image range to true min-max range")
    optgroup.add_option("--image_range_min",type="float",dest="image_range_min",
                      default=0.00, help="desired output intensity minimum") 
    optgroup.add_option("--image_range_max",type="float",dest="image_range_max",
                      default=0.00, help="desired output intensity maximum")
    optgroup.add_option("--volumeType",type="string",dest="vType",
                      default=None, help="datatype of mnc2 file")
    parser.add_option_group(optgroup)
    #sequence specific flag must be second command-line argument to specify sequence specific options and recon
    if (recontypestring is None):
        recontypestring = sys.argv[1]
    try:
        seqmodname=os.path.basename(recontypestring)
        seqdir=[[os.path.dirname(recontypestring)],None][os.path.dirname(recontypestring)=='']
        if seqdir is None:
            seqmod_specs = importlib.util.find_spec('mri_python.'+seqmodname)
        else:
            loader_details = (importlib.machinery.ExtensionFileLoader,importlib.machinery.EXTENSION_SUFFIXES)
            modfinder = importlib.machinery.FileFinder(seqdir, loader_details)
            seqmod_specs = modfinder.find_spec(seqmodname)

        if (isinstance(seqmod_specs.loader,zipimport.zipimporter)):
            seqmodule = seqmod_specs.loader.load_module('mri_python.'+seqmodname)
        else :
            seqmodule = importlib.util.module_from_spec(seqmod_specs)
            seqmod_specs.loader.exec_module(seqmodule)
    except ImportError:
        print("Could not import %s recon module"%recontypestring)
        seqmod_specs = importlib.util.find_spec('mri_python.'+DEFAULT_RECON_TYPE)
        seqmodule = importlib.util.module_from_spec(seqmod_specs)
        seqmod_specs.loader.exec_module(seqmodule)
        #f,filename,description = imp.find_module(DEFAULT_RECON_TYPE)
        #seqmodule = imp.load_module('seqmodule',f,filename,description)
    seqmodule.seq_specific_options(parser)
    #return options parser and seqmodule
    return parser,seqmodule


#----------------------------------------------------------------------
# top level program

if __name__ == '__main__':

    parser,seqmodule = generate_option_parser_and_seq_module()
    options, args = parser.parse_args()

    outputfile = args[-1]
    inputdirectory = args[-2]

    try:
        if not options.clobber and glob.glob(outputfile[:-4]+'*'): #os.path.exists(outputfile):
            raise FatalError("The --clobber option is needed to overwrite an existing file.")

        if not os.path.exists(inputdirectory):
            raise FatalError("Cannot locate input directory.")

    except FatalError as e:
        print('Error(%s):' % program_name, e.msg)
        raise SystemExit

    cmdstr = "renice 15 %d" % os.getpid()
    os.system(cmdstr)

    #determine if it is a bruker or varian directory and file format
    #then load acquisition details
    if (glob.glob(os.path.join(inputdirectory,options.procpar_file_name))):
        inputAcq = vrf.VarianAcquisition(inputdirectory,
                                         procpar_file=options.procpar_file_name)
    elif (glob.glob(os.path.join(inputdirectory,options.method_file_name)) and 
          glob.glob(os.path.join(inputdirectory,options.acqp_file_name))):
        inputAcq = brf.BrukerAcquisition(inputdirectory,
                                         method_name=options.method_file_name,
                                         acqp_name=options.acqp_file_name,
                                         fid_name=options.fid_file_name)
    else:
        print('Error: did not recognize input directory as Varian or Bruker format.')
        raise SystemExit

    #initialize seq specific data
    seqrec = seqmodule.seq_reconstruction(inputAcq,options,outputfile)

    #construct list of mice to reconstruct
    if (options.mouse_list):
        mouse_list = [int(x) for x in options.mouse_list.split(',') if N.in1d(int(x),inputAcq.rcvrmouse_mapping)]
    else:
        mouse_list = N.unique(inputAcq.rcvrmouse_mapping)  #range(inputAcq.nmice)

    #perform any global calibration or other analysis for all coils here
    if hasattr(seqrec,'pre_recon_processing'):
        seqrec.pre_recon_processing()

    #now perform recon, one mouse at a time
    for imouse in mouse_list:
        #get kspace
        seqrec.gen_kspace(imouse=imouse)
        #perform FFT
        if hasattr(seqrec,'recon'):
            seqrec.recon()
        else:
            default_recon(seqrec)
        #name file and output
        if (inputAcq.nmice>1):
            curroutputfile = outputfile[:-3]+str(imouse+1)+'.mnc'
        else:
            curroutputfile = outputfile
        if hasattr(seqrec,'output_to_file'):
            seqrec.output_to_file(curroutputfile,imouse)
        else:
            write_to_mnc_file(curroutputfile,seqrec.image_data,inputAcq,options,imouse=imouse)            
        #tidying
        del seqrec.kspace
        del seqrec.image_data
            


##for debugging in interactive python, it is convenient to define the command line options in a
##dummy class to feed into functions
class dummyopt():
    def __init__(self,petable=None,petable_ordered_pairs=False,fov_shift_ro="0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
                 fov_shift_pe1="0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",fov_shift_pe2="0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
                 outputreps=False,complexavg=False,petable_pe1=False,petable_pe2=False,fermi_ellipse=False,fermi_pecirc=False,
                 noshift_ppe=False,noshift_ppe2=False,large_data_recon=False,mouse_list=None,
                 nofft=False,fft1d=False,fft2d=False,fft3d=True,real=False,imag=False,phase=False,
                 image_range_min=0.0,image_range_max=0.0,max_range=False,vType=None,apowidth=0.8,
                 echoamp_alpha=0.02,phasecorr_data=None,phasecorr_table=None,phasecorr_plot=False,
                 no_echo_roshift=False,no_pe_phaseshift=False,no_Echo_shift_apply=False,no_echo_phaseshift=False,echo_ampl_only=False,
                 separate_multi_acq=False,fid_file_name="fid",alphamap=None,T1fit=False):
        self.petable=petable
        self.petable_ordered_pairs=petable_ordered_pairs
        self.fov_shift_ro=fov_shift_ro
        self.fov_shift_pe1=fov_shift_pe1
        self.fov_shift_pe2=fov_shift_pe2
        self.outputreps=outputreps
        self.complexavg=complexavg
        self.petable_pe1=petable_pe1
        self.petable_pe2=petable_pe2
        self.fermi_ellipse=fermi_ellipse
        self.fermi_pecirc=fermi_pecirc
        self.noshift_ppe=noshift_ppe
        self.noshift_ppe2=noshift_ppe2
        self.large_data_recon=large_data_recon
        self.mouse_list=mouse_list
        self.nofft=nofft
        self.fft1d=fft1d
        self.fft2d=fft2d
        self.fft3d=fft3d
        self.real=real
        self.imag=imag
        self.phase=phase
        self.image_range_min=image_range_min
        self.image_range_max=image_range_max
        self.max_range=max_range
        self.vType=vType
        self.apowidth=apowidth
        self.echoamp_alpha=echoamp_alpha
        self.phasecorr_data=phasecorr_data
        self.phasecorr_table=phasecorr_table
        self.phasecorr_plot=phasecorr_plot
        self.no_echo_roshift=no_echo_roshift
        self.no_echo_phaseshift=no_echo_phaseshift
        self.echo_ampl_only=echo_ampl_only
        self.no_pe_phaseshift=no_pe_phaseshift
        self.no_Echo_shift_apply=no_Echo_shift_apply
        self.separate_multi_acq=separate_multi_acq
        self.fid_file_name=fid_file_name
        self.alphamap=alphamap
        self.T1fit=T1fit

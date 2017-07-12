#!/usr/bin/env python
#
# cyltable_gen.py
#
# a script to generate petables for varian
# acquisition and reconstruction
#
# Created May 24, 2013

from sys import argv
import string
import os
import getopt
import re
from numpy import *
from numpy.random import shuffle
from optparse import OptionParser, Option, OptionValueError

program_name = 'cyltable_gen.py'

#----------------------------------------------------------------------------
# define program specific exception
class FatalError(Exception):
    def __init__(self,args=None):
        self.msg = args
#----------------------------------------------------------------------------

def cyl_acq_mask(nv2,nv,grappafov=1,tabledisc=1,etl=None,checkerboard=0):
    cartesian_distsqr = ((arange(0,nv2,grappafov)-nv2//2)/float(nv2//2))[:,newaxis]**2 + \
                        ((arange(0,nv,grappafov)-nv//2)/float(nv//2))[newaxis,:]**2
    if (checkerboard>0):
        allinds=meshgrid(arange(0,nv,grappafov)-nv//2,arange(0,nv2,grappafov)-nv2//2)
        gridtest = (allinds[0]%(2*grappafov) + allinds[1]%(2*grappafov))%(2*grappafov)
        moduloval = [0,grappafov][checkerboard-1]
        retval=nv*nv2 #just a number bigger than any of the others
        cartesian_distsqr = where(gridtest==moduloval,cartesian_distsqr,retval)
        ntotalpts=int( pi*nv*nv2/(4.0*2*grappafov**2) )
    else:
        ntotalpts=int( pi*nv*nv2/(4.0*grappafov**2) )
    if (ntotalpts>tabledisc):
        ntotalpts=tabledisc*( (ntotalpts//tabledisc) + [0,1][(ntotalpts%tabledisc)>0] )
    if not (etl is None):
        ntotalpts=etl*(ntotalpts/etl)
    cartesian_argsort = argsort(cartesian_distsqr.flat)[0:ntotalpts]
    acqmask = zeros((nv2,nv),bool)
    acqmask[grappafov*(cartesian_argsort//(nv//grappafov)),grappafov*(cartesian_argsort%(nv//grappafov))] = 1
    return acqmask

def pseudocyl_acq_mask(nv2,nv,tabledisc=1,etl=6,checkerboard=0):
    acqmask = cyl_acq_mask(nv2,nv,grappafov=1,tabledisc=tabledisc,checkerboard=checkerboard)
    ntotalpts = sum(acqmask)
    epts = 0
    stepval = [1,2][checkerboard>0]
    for j in range(nv2):
        npts_j = sum(acqmask[j,:])
        checkoffset = [0,j%2,(j+1)%2][checkerboard]
        if (npts_j<etl):
            acqmask[j,nv/2-stepval*(etl/2)+checkoffset:nv/2+stepval*(etl/2)+etl%2:stepval] = 1
            epts += etl-npts_j
        elif ((npts_j%etl)>0) and (epts>=0):
            q = etl*(npts_j/etl)
            acqmask[j,:] = 0
            acqmask[j,nv/2-stepval*(q/2)+checkoffset:nv/2+stepval*(q/2)+q%2:stepval] = 1
            epts += q-npts_j
        elif ((npts_j%etl)>0) and (epts<0):
            q = etl*(1+npts_j/etl)
            acqmask[j,:] = 0
            acqmask[j,nv/2-stepval*(q/2)+checkoffset:nv/2+stepval*(q/2)+q%2:stepval] = 1
            epts += q-npts_j
    return acqmask

def cart_acq_mask(nv2,nv,grappafov=1):
    acqmask = zeros((nv2,nv),bool)
    acqmask[::grappafov,::grappafov] = 1
    return acqmask

def grappa_acq_mask(nv2,nv,grappape2=16,grappape=16,grappa4prof=False,grappafov=8):
    acqmask = zeros((nv2,nv),bool)
    if (grappa4prof):
        acqmask[nv2/2-grappafov,(nv-grappape)/2-grappafov:(nv+grappape)/2-grappafov] = 1
        acqmask[(nv2-grappape2)/2-grappafov:(nv2+grappape2)/2-grappafov,nv/2-grappafov] = 1
        ind = arange(min([grappape,grappape2]))-min([grappape,grappape2])/2
        acqmask[nv2/2-grappafov+ind,nv/2-grappafov+ind] = 1
        acqmask[nv2/2-grappafov-ind,nv/2-grappafov+ind] = 1
    else:
        acqmask[(nv2-grappape2)/2:(nv2+grappape2)/2,(nv-grappape)/2:(nv+grappape)/2] = 1
    return acqmask

def distordwalk(acqmask):
    cpe1=0
    cpe2=0
    dpe1=-1
    dpe2=0
    nv=acqmask.shape[1]
    nv2=acqmask.shape[0]
    sorti1=empty(sum(acqmask),int)
    sorti2=empty(sum(acqmask),int)
    cpt=0
    for j in range(max(nv,nv2)**2):
        if (-nv/2<cpe1<=nv/2) and (-nv2/2<cpe2<=nv2/2):
            if (acqmask[nv2/2+cpe2,nv/2+cpe1]):
                sorti1[cpt]=cpe1
                sorti2[cpt]=cpe2
                cpt+=1
        cpe1,cpe2 = cpe1+dpe1,cpe2+dpe2
        if (cpe1==cpe2) or (cpe1<0 and cpe1==(-cpe2-1)) or (cpe1>0 and cpe1==-cpe2):
            dpe1,dpe2=-dpe2,dpe1
    return sorti2,sorti1

def angwalk(acqmask,etl=1,interleave=True,angdist=True,interleaveorder=None,nv=None,nv2=None):
    i1,i2 = nonzero(acqmask)
    angvals = arctan2( (i1-acqmask.shape[1]/2+1) , (i2-acqmask.shape[0]/2+1) )
    dist = sqrt( (i1-acqmask.shape[1]/2+1)**2 + (i2-acqmask.shape[0]/2+1)**2 )
    inds = argsort(dist)
    dist = dist[inds]
    angvals = angvals[inds]
    i1 = i1[inds]; i2=i2[inds]
    for j in range(etl):
        inds = argsort(angvals[j*len(dist)/etl:(j+1)*len(dist)/etl])
        i1[j*len(dist)/etl:(j+1)*len(dist)/etl] = i1[inds+j*len(dist)/etl]
        i2[j*len(dist)/etl:(j+1)*len(dist)/etl] = i2[inds+j*len(dist)/etl]
    if (angdist): #re-sort within each etl by pe2 distance (but keep angles matched)
        inds = argsort(i2[0:len(dist)/etl])
        for j in range(etl):
            i1[j*len(dist)/etl:(j+1)*len(dist)/etl] = i1[j*len(dist)/etl:(j+1)*len(dist)/etl][inds]
            i2[j*len(dist)/etl:(j+1)*len(dist)/etl] = i2[j*len(dist)/etl:(j+1)*len(dist)/etl][inds]
    if (interleave):
        sorti1 = empty(i1.shape,int)
        sorti2 = empty(i2.shape,int)
        if (interleaveorder==None):
            interleaveorder=arange(etl)
        for j in range(len(i1)/etl):
            sorti1[j*etl:(j+1)*etl] = i1[j::len(i1)/etl][interleaveorder]
            sorti2[j*etl:(j+1)*etl] = i2[j::len(i2)/etl][interleaveorder]
    else:
        sorti1=i1
        sorti2=i2
    if (nv==None): nv=acqmask.shape[1]
    if (nv2==None): nv2=acqmask.shape[0]
    sorti1 -= nv/2-1
    sorti2 -= nv2/2-1
    return sorti2,sorti1

def cylaxisdistwalk(acqmask,etl=1,k0echo=0,nv=None,nv2=None,axisangle=0):
    i1,i2 = nonzero(acqmask)
    dist = sqrt( (i1-acqmask.shape[1]/2+1)**2 + (i2-acqmask.shape[0]/2+1)**2 )
    #sort first by i1'
    inds = argsort(i1*cos(pi*(axisangle)/180.0) - \
                   i2*sin(pi*(axisangle)/180.0))
    dist = dist[inds]
    i1 = i1[inds]; i2=i2[inds]
    #sort second by dist on each half i1'
    for j in range(2):
        inds = argsort( dist[j*len(dist)/2:(j+1)*len(dist)/2] )
        if (j==0):
            inds = inds[::-1]
        dist[j*len(dist)/2:(j+1)*len(dist)/2] = dist[inds+j*len(dist)/2]
        i1[j*len(dist)/2:(j+1)*len(dist)/2] = i1[inds+j*len(dist)/2]
        i2[j*len(dist)/2:(j+1)*len(dist)/2] = i2[inds+j*len(dist)/2]
    #roll into annuli
    if (etl%2): #odd number of echoes
        i1 = roll(i1,(len(i1)/etl)*(k0echo-etl/2))
        i2 = roll(i2,(len(i2)/etl)*(k0echo-etl/2))
    else:       #even number of echoes
        i1 = roll(i1,(len(i1)/etl)*(k0echo-etl/2+1)-len(i1)/(2*etl))
        i2 = roll(i2,(len(i2)/etl)*(k0echo-etl/2+1)-len(i2)/(2*etl))
    #sort last by i2' within annuli
    for j in range(etl):
        inds = argsort( i1[j*len(dist)/etl:(j+1)*len(dist)/etl]*sin(pi*(axisangle)/180.0) + \
                        i2[j*len(dist)/etl:(j+1)*len(dist)/etl]*cos(pi*(axisangle)/180.0) )
        i1[j*len(dist)/etl:(j+1)*len(dist)/etl] = i1[inds+j*len(dist)/etl]
        i2[j*len(dist)/etl:(j+1)*len(dist)/etl] = i2[inds+j*len(dist)/etl]
    #reorder based on specified interleaveorder
    interleaveorder = arange(etl)
    sorti1 = empty(i1.shape,int)
    sorti2 = empty(i2.shape,int)
    for j in range(len(i1)/etl):
        sorti1[j*etl:(j+1)*etl] = i1[j::len(i1)/etl][interleaveorder]
        sorti2[j*etl:(j+1)*etl] = i2[j::len(i2)/etl][interleaveorder]
    if (nv==None): nv=acqmask.shape[1]
    if (nv2==None): nv2=acqmask.shape[0]
    sorti1 -= nv/2-1
    sorti2 -= nv2/2-1
    return sorti2,sorti1

def pseudocyl(acqmask,etl=1,k0echo=0,nv=None,nv2=None):
    npts = len(nonzero(acqmask)[0])
    i2 = empty((npts,),int)
    i1 = empty((npts,),int)
    if (nv==None): nv=acqmask.shape[1]
    if (nv2==None): nv2=acqmask.shape[0]
    cpt = 0
    for j in range(acqmask.shape[0]):
        inds = nonzero(acqmask[j,:])[0]
        gsize = len(inds)/etl
        if (gsize==0): continue
        if ((etl%2)==0):
            egrp = roll(arange(len(inds))/gsize,(etl/2-k0echo)*gsize+gsize/2)
        else:
            egrp = roll(arange(len(inds))/gsize,(etl/2-k0echo)*gsize)
        i_eord = empty((len(inds)/etl,etl),int)
        for k in range(etl):
            i_eord[:,k]=sort(inds[egrp==k])
        i_eord = ravel(i_eord)
        i1[cpt:cpt+len(i_eord)] = i_eord.copy()-acqmask.shape[1]/2+1
        i2[cpt:cpt+len(i_eord)] = j-acqmask.shape[0]/2+1
        cpt += len(i_eord)
    return i2,i1

def cyl_concentr_ellipses(acqmask,etl=1,k0echo=None,nv=None,nv2=None):
    npts = len(nonzero(acqmask)[0])
    NTR = npts/etl
    if (nv==None): nv=acqmask.shape[1]
    if (nv2==None): nv2=acqmask.shape[0]
    if (k0echo is None):
        k0echo=etl/2
    i0,i1=nonzero(acqmask)
    i0=i0-acqmask.shape[0]/2
    i1=i1-acqmask.shape[1]/2
    etl_assign=zeros((nv2,nv),int)
    #middle ROI
    E_val=i0**2/float(nv2/2)**2+i1**2/(nv/float(2*etl))**2
    cinds=argsort(E_val)[0:NTR]
    etl_assign[i0[cinds]+nv2/2,i1[cinds]+nv/2]=etl/2
    #symmetric rings
    Neginds=nonzero(i1<0)[0]
    Posinds=nonzero(i1>=0)[0]
    j=1
    while (etl/2+j<etl):
        E_val=i0**2/float(nv2/2)**2+i1**2/((j+0.5)*nv/float(2*etl))**2
        E_val[etl_assign[i0+nv2/2,i1+nv/2]>0]=nv2**2+nv**2 #artificial assignment to already "acquired" points
        #negative side
        cinds=Neginds[argsort(E_val[Neginds])[0:NTR]]
        etl_assign[i0[cinds]+nv2/2,i1[cinds]+nv/2]=etl/2-j
        #positive side
        cinds=Posinds[argsort(E_val[Posinds])[0:NTR]]
        etl_assign[i0[cinds]+nv2/2,i1[cinds]+nv/2]=etl/2+j
        j=j+1
    #for even etl, one last ring that is split on both sides
    if (etl%2==0):
        etl_assign=where((etl_assign==0)&(acqmask>0),etl,etl_assign)
    #perform the ordering
    i0_ordered=zeros(len(i0),int)
    i1_ordered=zeros(len(i1),int)
    for j in range(etl):
        ci0,ci1=nonzero(etl_assign==(j+1))
        roi_ord=argsort(ci0)
        i0_ordered[j::etl]=ci0[roi_ord]-nv2/2
        i1_ordered[j::etl]=ci1[roi_ord]-nv/2
    #shuffle for desired centre
    for j in range(NTR):
        i0_ordered[j*etl:(j+1)*etl]=roll(i0_ordered[j*etl:(j+1)*etl],(k0echo-etl/2))
        i1_ordered[j*etl:(j+1)*etl]=roll(i1_ordered[j*etl:(j+1)*etl],(k0echo-etl/2))
    return i0_ordered,i1_ordered

def mp2rage_replicate(i2,i1,etl,mp2rage_echoshare):
    i2mp2=zeros(2*(len(i2)-mp2rage_echoshare*etl),int)
    i1mp2=zeros(2*(len(i1)-mp2rage_echoshare*etl),int)
    etlnew=2*(etl-mp2rage_echoshare)
    for j in range(len(i2)/etl):
        i2mp2[j*etlnew:(j*etlnew+etl-mp2rage_echoshare)]=i2[(j*etl+mp2rage_echoshare):(j+1)*etl]
        i2mp2[(j*etlnew+etl-mp2rage_echoshare):(j*etlnew+etlnew)]=i2[j*etl:((j+1)*etl-mp2rage_echoshare)]
        i1mp2[j*etlnew:(j*etlnew+etl-mp2rage_echoshare)]=i1[(j*etl+mp2rage_echoshare):(j+1)*etl]
        i1mp2[(j*etlnew+etl-mp2rage_echoshare):(j*etlnew+etlnew)]=i1[j*etl:((j+1)*etl-mp2rage_echoshare)]        
    return i2mp2,i1mp2

def feather45(i2,i1,etl):
    i12_flag = (((i1+i2)%2)==1)
    ecurr = (arange(len(i1))%etl)
    etarg = (ecurr+1)%etl
    for j in range(len(i12_flag)):
        if (not i12_flag[j]): continue
        etest = (etarg[(j+1):]==ecurr[j])*(i12_flag[(j+1):])
        if not etest.any():
            etest = ((etarg[(j+1):]%2)==(ecurr[j]%2))*(i12_flag[(j+1):])     
        if etest.any(): 
            cind = j+1+nonzero(etest)[0][0]
            temp = i1[j]; i1[j]=i1[cind]; i1[cind]=temp
            temp = i2[j]; i2[j]=i2[cind]; i2[cind]=temp
            temp = etarg[j]; etarg[j]=etarg[cind]; etarg[cind]=temp            
            i12_flag[j] = False
    return i2,i1

def NTRfeather(i2,i1,NTRfeather):
    i1mod = i1.copy()
    i2mod = i2.copy()
    for j in range(len(i2)/(2*NTRfeather)):
        inds1 = i1[j*2*NTRfeather:(j+1)*2*NTRfeather]
        i1mod[j*2*NTRfeather:j*2*NTRfeather+NTRfeather] = inds1[0:2*NTRfeather:2]
        i1mod[j*2*NTRfeather+NTRfeather:(j+1)*2*NTRfeather] = inds1[1:2*NTRfeather:2]
        inds2 = i2[j*2*NTRfeather:(j+1)*2*NTRfeather]
        i2mod[j*2*NTRfeather:j*2*NTRfeather+NTRfeather] = inds2[0:2*NTRfeather:2]
        i2mod[j*2*NTRfeather+NTRfeather:(j+1)*2*NTRfeather] = inds2[1:2*NTRfeather:2]
    return i2mod,i1mod

#def gen_phase_corr(i1,i2,TRundersample,etl):
#    I_step = etl*TRundersample
#    #pc_i1 = zeros([1+(etl-1)*(len(i1)/I_step),etl],int)
#    #pc_i2 = zeros([1+(etl-1)*(len(i2)/I_step),etl],int)
#    pc_i1 = zeros([1+(etl)*(len(i1)/I_step),etl],int)
#    pc_i2 = zeros([1+(etl)*(len(i2)/I_step),etl],int)
#    for cTR in range(0,I_step*(len(i1)/I_step),I_step):
#        c_pcind = 1+(etl)*(cTR/I_step)
#        for j in range(etl):
#            pc_i1[c_pcind+j,0:(j+1)] = i1[cTR:cTR+j+1]
#            pc_i2[c_pcind+j,0:(j+1)] = i2[cTR:cTR+j+1]
#    pc_i1 = reshape(pc_i1,(pc_i1.shape[0]*pc_i1.shape[1],))
#    pc_i2 = reshape(pc_i2,(pc_i2.shape[0]*pc_i2.shape[1],))
#    return pc_i1,pc_i2

def gen_phase_corr(matrix_pe1,matrix_pe2,matundersample,etl):
    I_step = matundersample
    Npe1 = matrix_pe1/I_step
    Npe2 = matrix_pe2/I_step
    Npe1pe2 = min([Npe1,Npe2])
    pc_i1 = zeros([Npe1+Npe2+2*Npe1pe2,etl],int)
    pc_i2 = zeros([Npe1+Npe2+2*Npe1pe2,etl],int)
    #pe1 range, no pe 2
    pc_i1[0:Npe1,0] = arange(0,Npe1*I_step,I_step)-Npe1*I_step/2
    #no pe1, pe2 range
    pc_i2[Npe1:Npe1+Npe2,0] = arange(0,Npe2*I_step,I_step)-Npe2*I_step/2
    #pe1 & pe2 range
    pc_i1[Npe1+Npe2:Npe1+Npe2+Npe1pe2,0] = arange(0,Npe1pe2*I_step,I_step)-Npe1pe2*I_step/2
    pc_i2[Npe1+Npe2:Npe1+Npe2+Npe1pe2,0] = arange(0,Npe1pe2*I_step,I_step)-Npe1pe2*I_step/2
    #pe1 & -pe2 range
    pc_i1[Npe1+Npe2+Npe1pe2:Npe1+Npe2+2*Npe1pe2,0] = arange(0,Npe1pe2*I_step,I_step)-Npe1pe2*I_step/2
    pc_i2[Npe1+Npe2+Npe1pe2:Npe1+Npe2+2*Npe1pe2,0] = -1*(arange(0,Npe1pe2*I_step,I_step)-Npe1pe2*I_step/2)
    #flatten for output
    pc_i1 = reshape(pc_i1,(pc_i1.shape[0]*pc_i1.shape[1],))
    pc_i2 = reshape(pc_i2,(pc_i2.shape[0]*pc_i2.shape[1],))
    return pc_i1,pc_i2

def gen_phase_corrfull(matrix_pe1,matrix_pe2,matundersample,etl):
    I_step = matundersample
    Npe1 = matrix_pe1/I_step
    Npe2 = matrix_pe2/I_step
    Npe1pe2 = min([Npe1,Npe2])
    pc_i1 = zeros([Npe1*Npe2,etl],int)
    pc_i2 = zeros([Npe1*Npe2,etl],int)
    for j in range(Npe2):
        pc_i1[j*Npe1:(j+1)*Npe1,0]=arange(0,Npe1*I_step,I_step)-Npe1*I_step/2      
        pc_i2[j*Npe1:(j+1)*Npe1,0]=j*I_step-Npe2*I_step/2
    #flatten for output
    pc_i1 = reshape(pc_i1,(pc_i1.shape[0]*pc_i1.shape[1],))
    pc_i2 = reshape(pc_i2,(pc_i2.shape[0]*pc_i2.shape[1],))
    return pc_i1,pc_i2
    
def petable_file_output(i1,i2,nf,ntables,outputfile,linelim=8,appendnfni=True):
    file_num_format = "_%d" 
    ni = len(i1)/(nf*ntables)
    ptspertable = len(i1)/ntables
    if (appendnfni):
        outputfile = outputfile + "_nf%d_ni%d"%(nf,ni)
    full_fh = open(outputfile,'w')
    full_fh.write("t1 = \n")
    file_fh_list=[]
    multifid = (ntables>1)
    for j in range(ntables):
        if (multifid):
            curr_fh = open(outputfile+file_num_format%j,'w')
            curr_fh.write("t1 = \n")
            file_fh_list.append(curr_fh)
        for k in range(ni):
            for q in range(nf):
                full_fh.write("    %d"%i1[j*ptspertable+k*nf+q])
                if (multifid): curr_fh.write("    %d"%i1[j*ptspertable+k*nf+q])
                if ((q+1)%linelim==0) and (q<nf-1):
                    full_fh.write("\n")
                    if (multifid): curr_fh.write("\n")
            full_fh.write("\n")
            if (multifid): curr_fh.write("\n")
    full_fh.write("t2 = \n")
    for j in range(ntables):
        if (multifid):
            curr_fh = file_fh_list[j]
            curr_fh.write("t2 = \n")
        for k in range(ni):
            for q in range(nf):
                full_fh.write("    %d"%i2[j*ptspertable+k*nf+q])
                if (multifid): curr_fh.write("    %d"%i2[j*ptspertable+k*nf+q])
                if ((q+1)%linelim==0) and (q<nf-1):
                    full_fh.write("\n")
                    if (multifid): curr_fh.write("\n")
            full_fh.write("\n")
            if (multifid): curr_fh.write("\n")
        if (multifid): curr_fh.close()
    full_fh.close()
    return None

#----------------------------------------------------------------------
# top level program

if __name__ == '__main__':

    usage = """%s <matrix1> <matrix2> <output file>
   or  %s --help
   
%s is a script for generating petable files for varian cylindrical ge image acquisition
and reconstruction.
"""
    usage = usage % ((program_name, )*3)

    parser = OptionParser(usage)
    parser.add_option("--clobber", action="store_true", dest="clobber",
                       default=0, help="overwrite output file")
    parser.add_option("--bruker",action="store_true",dest="bruker",
                       default=0, help="table for bruker acquistion (no need for multiple fid files or ni*nf descretization)")
    parser.add_option("--cylindrical", action="store_true", dest="cylindrical",
                       default=0, help="shave corners of k-space in PE1 and PE2")
    parser.add_option("--distorder", action="store_true", dest="distorder", default=0,
                      help="order elements from low to high (square spiral ordering)")
    parser.add_option("--angorder", action="store_true", dest="angorder", default=0,
                      help="order elements by angle in PE2-PE1 plane")
    parser.add_option("--angdist", action="store_true", dest="angdist", default=0,
                      help="order elements by angle in PE2-PE1 plane and distance in PE2")
    parser.add_option("--cylaxisdist", action="store_true", dest="cylaxisdist", default=0,
                      help="order elements so etl runs preferentially along an axis")
    parser.add_option("--cylaxisangle",type="float",dest="cylaxisangle",
                       default=0.0, help="preferred angle in degrees for 'fast' acquisition (default: 0.0 = PE1 axis)")
    parser.add_option("--checkerboard",type="int",dest="checkerboard",
                       default=0, help="grab only encodes on a checkerboard (0=False,1=black squares, 2=white squares")    
    parser.add_option("--pseudocyl",action="store_true",dest="pseudocyl",default=0,
                      help="linearly arrange etl along pe1 axis within cylinder, similar to cartesian images")
    parser.add_option("--pseudocyl_pe2",action="store_true",dest="pseudocyl_pe2",default=0,
                      help="linearly arrange etl along pe2 axis within cylinder, similar to cartesian images")
    parser.add_option("--cyl_concentric_ellipses",action="store_true",dest="cyl_concentr_ellipses",default=0,
                      help="sort multiple 'echoes' by concentric ellipses (major axis on pe2) with etl divisions")
    parser.add_option("--randomize", action="store_true", dest="randomize", default=0,
                      help="randomize the order of acquisition")
    parser.add_option("--nf",type="int",dest="nf",
                       default=1, help="desired nf for each computed table")
    parser.add_option("--grappafov",type="int",dest="grappafov",
                       default=1, help="controls size of full FOV relative multi-mouse FOV")
    parser.add_option("--grappape",type="int",dest="grappape",
                       default=0, help="PE for grappa data")
    parser.add_option("--grappape2",type="int",dest="grappape2",
                       default=0, help="PE2 for grappa data")
    parser.add_option("--grappa4prof",action="store_true",dest="grappa4prof",
                       default=0, help="Acquire only 4 profiles for grappa data instead of full k-space")
    parser.add_option("--individ_table_limit",type="int",dest="individ_table_limit",
                       default=4096, help="limit of number of elements in each table")
    parser.add_option("--nreps",type="int",dest="nreps",
                       default=1, help="repeat a single table specified number of times")
    parser.add_option("--repk0",type="int",dest="repk0",
                       default=-1, help="repeat k0 acquisition approximately every # of TRs")
    parser.add_option("--etl",type="int",dest="etl",
                       default=1, help="Echo train length")
    parser.add_option("--etlorder",type="string",dest="etlorder",
                       default=None, help="Echo train kspace order (0..etl-1)")
    parser.add_option("--k0echo",type="int",dest="k0echo",
                       default=1, help="echo number for centre of k-space in cylpe2dist option")
    parser.add_option("--feather45", action="store_true", dest="feather45",
                       default=0, help="feather echoes at 45 degree angle to trade-off ripple vs ghost")
    parser.add_option("--NTRfeather",type="int",dest="NTRfeather",
                       default=0, help="interleave every N phase encode points before outputing final table (intended for use with gradient-echo sequences)")
    parser.add_option("--output_phase_corr_file", action="store_true", dest="output_phase_corr_file",
                       default=0, help="output a reduced sampling phase correction file (covers principle axes and diagonals only)")
    parser.add_option("--output_phase_corr_filefull", action="store_true", dest="output_phase_corr_filefull",
                       default=0, help="output a phase correction file that covers all of 2D cartesian space")
    parser.add_option("--phase_corr_sampling",type="int",dest="phase_corr_sampling",
                       default=1, help="Undersampling for phase correction file")
    parser.add_option("--mp2rage_rep", action="store_true", dest="mp2rage_rep",
                       default=0, help="double etl to obtain table suitable for mp2rage_rep")
    parser.add_option("--mp2rage_echoshare",type="int",dest="mp2rage_echoshare",
                       default=0, help="number of shared echos in mp2rage train (default: 0)")

    options, args = parser.parse_args()

    outputfile = args[-1]
    try:
        if not options.clobber and os.path.exists(outputfile):
            raise FatalError("The --clobber option is needed to overwrite an existing file.")
        if (options.nf%options.etl!=0):
            raise FatalError("nf and etl are not compatible.")
        matrix_pe1 = int(args[-3])
        matrix_pe2 = int(args[-2])
    except FatalError as e:
        print('Error(%s):' % program_name, e.msg)
        raise SystemExit

    #get etl order if specified
    if ((options.etlorder) and (options.etl>1)):
        etl_order_array = argsort( array(options.etlorder.split(","),int) )
    else:
        etl_order_array = arange(options.etl)

    #set the acquisition mask
    if (options.bruker):
        w_tablimit = options.nf*(512**2//options.nf) #bruker table limit is 1e6?? set it to 512**2 here
    else:
        w_tablimit = options.nf*(options.individ_table_limit//options.nf)
    if ((options.pseudocyl) or (options.pseudocyl_pe2)):
        if (options.pseudocyl):
            acqmask = pseudocyl_acq_mask(matrix_pe2,
                                         matrix_pe1,
                                         tabledisc=w_tablimit,
                                         etl=options.etl,
                                         checkerboard=options.checkerboard)
        else: #pseudocyl_pe2, same as pseudocyl but with axes swapped
            acqmask = transpose( pseudocyl_acq_mask(matrix_pe1,
                                                    matrix_pe2,
                                                    tabledisc=w_tablimit,
                                                    etl=options.etl,
                                                    checkerboard=options.checkerboard) )
    elif (options.cylindrical):
            acqmask = cyl_acq_mask(matrix_pe2*options.grappafov,
                                   matrix_pe1*options.grappafov,
                                   grappafov=options.grappafov,
                                   tabledisc=w_tablimit,
                                   etl=options.etl,
                                   checkerboard=options.checkerboard)
    else:
        acqmask = cart_acq_mask(matrix_pe2*options.grappafov,
                                matrix_pe1*options.grappafov,
                                grappafov=options.grappafov)

    if (options.grappafov>1):
        grappamask = grappa_acq_mask(matrix_pe2*options.grappafov,
                                     matrix_pe1*options.grappafov,
                                     grappape2=options.grappape2,
                                     grappape=options.grappape,
                                     grappa4prof=options.grappa4prof,
                                     grappafov=options.grappafov)
        acqmask = acqmask + grappamask

    
    #sort the array into an index list based on the desired order
    if (options.distorder):
        i2,i1 = distordwalk(acqmask)
    elif ((options.angorder) or (options.angdist)):
        i2,i1 = angwalk(acqmask,etl=options.etl,angdist=options.angdist,interleaveorder=etl_order_array,nv=matrix_pe1,nv2=matrix_pe2)
    elif (options.cylaxisdist):
        k0echo = options.k0echo-1
        i2,i1 = cylaxisdistwalk(acqmask,etl=options.etl,k0echo=options.k0echo,nv=matrix_pe1,nv2=matrix_pe2,axisangle=options.cylaxisangle)
    elif (options.pseudocyl):
        i2,i1 = pseudocyl(acqmask,etl=options.etl,k0echo=options.k0echo,nv=matrix_pe1,nv2=matrix_pe2)
    elif (options.pseudocyl_pe2): #same as pseudocyl but with axes swapped
        i1,i2 = pseudocyl(transpose(acqmask),etl=options.etl,k0echo=options.k0echo,nv=matrix_pe2,nv2=matrix_pe1)
    elif (options.cyl_concentr_ellipses):
        i2,i1 = cyl_concentr_ellipses(acqmask,etl=options.etl,k0echo=options.k0echo,nv=matrix_pe1,nv2=matrix_pe2)
    elif (options.randomize):
        i2,i1 = nonzero(acqmask)
        i2 = i2 - (matrix_pe2/2-1)*options.grappafov
        i1 = i1 - (matrix_pe1/2-1)*options.grappafov
        inds = arange(i2.shape[0])
        shuffle(inds)
        i2=i2[inds]
        i1=i1[inds]
    else: #default is a raster ordering
        i2,i1 = nonzero(acqmask)
        i2 = i2 - (matrix_pe2/2-1)*options.grappafov
        i1 = i1 - (matrix_pe1/2-1)*options.grappafov


    #if requested, stripe even/odd echoes on 45 degree angle
    if (options.feather45):
        i2,i1 = feather45(i2,i1,options.etl)

    if (options.NTRfeather>0):
        i2,i1 = NTRfeather(i2,i1,options.NTRfeather)

    if (options.mp2rage_rep):
        i2,i1 = mp2rage_replicate(i2,i1,options.etl,options.mp2rage_echoshare)

    #repeat the array in series after itself if requested
    if (options.nreps>1):
        npts = len(i1)
        for j in range(1,options.nreps):
            i1=append(i1,i1[0:npts])
            i2=append(i2,i2[0:npts])

    #insert repeated k0 acquisitions if requested
    if (options.repk0>0):
        for j in range(len(i2)/(options.repk0*options.etl)):
            i1=insert(i1,j*options.repk0*ones((options.etl,),int),zeros((options.etl,),int))
            i2=insert(i2,j*options.repk0*ones((options.etl,),int),zeros((options.etl,),int))

    if (not options.bruker):
        #now determine the adjusted number of points required
        nf = options.nf
        w_tablimit = nf*(options.individ_table_limit/nf)
        ntables = len(i2)/w_tablimit
        if (ntables==0): ntables = 1
        totalpts = nf*ntables*((len(i2)/(nf*ntables))+[0,1][len(i2)%(nf*ntables)>0])
        ni = totalpts/(nf*ntables)
        excesspts = totalpts-len(i2)
        #for excess pts, re-acquire (0,0) at approximately equal spacing throughout the scan
        #or tack on zeros at end if options.repk0 already specified
        #or just grab extra repeats of existing scan
        if ((excesspts>0) and (options.repk0<0)):
            ex_spacing = totalpts/excesspts
            for j in range(excesspts):
                i1=insert(i1,j*ex_spacing*ones((options.etl,),int),zeros((options.etl,),int))
                i2=insert(i2,j*ex_spacing*ones((options.etl,),int),zeros((options.etl,),int))
        elif (excesspts>0):
            i1=append(i1,zeros((excesspts,),int))
            i2=append(i2,zeros((excesspts,),int))
        ntables = len(i1)/(nf*ni)
    else:
        if ((len(i1)%options.nf)==0):
            nf = max([options.nf,options.etl])
        else:
            nf = options.etl
        ni = len(i1)/nf 
        ntables = 1

    #output to file
    petable_file_output(i1,i2,nf,ntables,outputfile,linelim=[8,options.etl][options.etl>1],appendnfni=(not options.bruker))

    #advise user of parameters needed for scan
    print("nf = %d"%options.nf)
    print("ni = %d"%ni)
    print("nfid = %d"%ntables)
    print("nTR = %d"%len(i2))
    print("kspace %% = %f"%(100*len(i2)/float(matrix_pe1*matrix_pe2)))

    #generate phase correction file if needed
    if ((options.output_phase_corr_file) or (options.output_phase_corr_filefull)):
        if (not options.output_phase_corr_filefull):
            pc_i1,pc_i2 = gen_phase_corr(matrix_pe1,matrix_pe2,options.phase_corr_sampling,options.etl)
        else:
            pc_i1,pc_i2 = gen_phase_corrfull(matrix_pe1,matrix_pe2,options.phase_corr_sampling,options.etl)
        nphasetables=1
        petable_file_output(pc_i1,pc_i2,options.nf,nphasetables,outputfile+'_phasecorr',linelim=options.etl,appendnfni=(not options.bruker))
        print("phasecorr nf = %d"%options.nf)
        print("phasecorr ni = %d"%(len(pc_i2)/options.nf))
        print("phasecorr nTR = %d"%(len(pc_i2)/options.etl))


#class options_struct:
#    def __init__(self):
#        self.clobber=False
#        self.cylindrical=True
#        self.distorder=True
#        self.randomize=False
#        self.nf=32
#        self.grappafov=8
#        self.grappape=32
#        self.grappape2=32
#        self.individ_table_limit=4096
#        self.nreps=2
# 
#options=options_struct()

#class options_struct:
#    def __init__(self):
#        self.clobber=False
#        self.cylindrical=True
#        self.distorder=False
#        self.randomize=False
#        self.cylpe2dist=True
#        self.etl=6
#        self.nf=60
#        self.grappafov=1
#        self.grappape=0
#        self.grappape2=0
#        self.individ_table_limit=2048
#        self.nreps=1
#        self.k0echo=3
#        self.feather45=True
# 
#options=options_struct()

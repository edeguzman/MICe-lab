#!/bin/bash 

module purge
module load CellProfiler
#module list
cellprofiler $@

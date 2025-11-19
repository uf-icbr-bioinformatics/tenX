#!/bin/bash

TENX_HOME=/apps/dibig_tools/tenX/

name=$1
if [[ -z $name ]]; then name=Seurat; fi

module load conda/23.10

mamba create -p ${TENX_HOME}/${name} -c conda-forge r-base -y
mamba activate ${TENX_HOME}/${name}
mamba install -c conda-forge r-essentials hdf5

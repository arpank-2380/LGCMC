#!/bin/bash
root_dir=`pwd`
binary_location=${root_dir}/bin
compiler_module='ic111_20100414'     # You may have to change this to look for your intel fortran compiler
module load ${compiler_module}
cd ${root_dir}/source
   make clean
   make all

if [[ -d ${binary_location} ]]; then
    mv lgcmc.x ${binary_location}
else
    mkdir ${binary_location}
    mv lgcmc.x ${binary_location}
fi

make clean
cd $root_dir

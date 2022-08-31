#! user/bin/bash 

module load netcdf
ml

export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH
export ECCIROOT=$(pwd)
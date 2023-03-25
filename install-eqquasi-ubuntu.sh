#! user/bin/bash 

module load netcdf
ml

export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
ln -s /usr/lib/x86_64-linux-gnu/libblas.so.3 /usr/lib/x86_64-linux-gnu/libblas.so
ln -s /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/liblapack.so

export TACC_NETCDF_INC=/usr/include
export TACC_NETCDF_LIB=/usr/lib/x86_64-linux-gnu
export FC=mpif90
export FLAG='-fopenmp -ffree-line-length-none'

cd src
make
cd ..
mkdir bin
mv src/eqquasi bin

export EQQUASIROOT=$(pwd)
export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH

chmod -R 755 scripts

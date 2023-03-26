#! user/bin/bash 

# The shell script is to set up environments for EQquasi and 
#	install it. It will call the makefile inside src/ and generate 
#	an executable eqquasi and move it to bin/.

# Currently, the machines supported are:
#	ls6:    Lonestar6 at TACC
# 	ubuntu: Ubuntu 22.04

echo "Users need to specify the ENV VAR MACHINE"

MACHINE=ubuntu # ls6/ubuntu

if [ $MACHINE == "ls6" ]; then 
	echo "Installing EQquasi on Lonestar6 at TACC ... ..."
	module load netcdf
	ml
	echo "NETCDF INC and LIB PATH"
	echo $TACC_NETCDF_INC
	echo $TACC_NETCDF_LIB
	
elif [ $MACHINE == "ubuntu" ]; then 
	echo "Installing EQquasi on Ubuntu 22.04 ... ..."
	export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
	ln -sf /usr/lib/x86_64-linux-gnu/libblas.so.3 /usr/lib/x86_64-linux-gnu/libblas.so
	ln -sf /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/liblapack.so
fi 

cd src
make
cd ..
mkdir bin
mv src/eqquasi bin

export EQQUASIROOT=$(pwd)
export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH

chmod -R 755 scripts

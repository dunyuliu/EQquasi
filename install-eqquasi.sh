#! user/bin/bash 

# The shell script is to set up environments for EQquasi and 
#	install it. It will call the makefile inside src/ and generate 
#	an executable eqquasi and move it to bin/.

# Currently, the machines supported are:
#	ls6:	Lonestar6 at TACC
#	 ubuntu: Ubuntu 22.04

# Usage: install-eqquasi [-h] [-m] Machine_name

while getopts "hm:" OPTION; do
    case $OPTION in
        m)
            MACH=$OPTARG
        h)
            echo "Full command: install-eqquasi [-h] [-m] Machine_name       "
            echo "Usage:                                                     "
            echo " source install-eqquasi -h                                 "
            echo "      Display this help message                            "
            echo " source install-eqquasi -m ls6                             "
            echo "      Install EQquasi on Lonestar6 at TACC                 "
            echo "                                                           "
            echo "Currently supported machines:                              "
            echo "  ls6/ubuntu                                               "
            ;;
    esac
done 
            ;;
# for arg in "$@"; do 
    # case $arg in 
        # --help)
            # echo "Usage: install-eqquasi [--help] MACH_name"
            # echo "Options:"
            # echo "  --help  Display this help message"
            # echo "Currently supported MACH name:"
            # echo "  ls6/ubuntu"
            # echo "Example:"
            # echo "  source install-eqquasi ls6"
            # ;;
    # esac
# done 

if [ -n "$MACH" ]; then 
    export MACHINE=$MACH
    if [ $MACHINE == "ls6" ]; then 
        echo "Installing EQquasi on Lonestar6 at TACC ... ..."
        module load netcdf mumps
        ml
        echo "NETCDF INC and LIB PATH"
        echo $TACC_NETCDF_INC
        echo $TACC_NETCDF_LIB
        echo $TACC_MUMPS_INC
        echo $TACC_MUMPS_LIB
        
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
fi


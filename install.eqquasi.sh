#! /bin/bash 

# The shell script is to set up environments for EQquasi and 
#	install it. It will call the makefile inside src/ and generate 
#	an executable eqquasi and move it to bin/.

# Currently, the machines supported are:
#	ls6:	Lonestar6 at TACC
#	ubuntu: Ubuntu 22.04

# Usage: install-eqquasi [-h] [-m Machine_name] [-c Machine_name]
set -- "$@"

while getopts "hm:c:" OPTION; do
    case $OPTION in
        m)
            MACH=$OPTARG
            ;;
        c)
            MACH=$OPTARG
            CONFIG="True"
            ;;
        h)
            echo "Usage: ./install-eqquasi.sh [-h] [-m Machine_name] [-c Machine_name] "
            echo "                                                                     "
            echo "Examples:                                                            "
            echo "                                                                     "
            echo "./install-eqquasi.sh -h                                              "
            echo " -----Display this help message                                      "
            echo "                                                                     "
            echo "./install-eqquasi.sh -m ls6                                          "
            echo " -----Install EQquasi on Lonestar6 at TACC                           "
            echo "                                                                     "
            echo "./install-eqquasi.sh -c ubuntu                                       "
            echo " -----Simply set up envs for EQquasi without installation            "
            echo " -----on ubuntu                                                      "
            echo "                                                                     "
            echo "source install-eqquasi.sh                                            "
            echo " -----Activate ENV VAR EQQUASIROOT and add exes to PATH              "
            echo "                                                                     "
            echo "Currently supported machines include:                                "
            echo " ls6/ubuntu                                                          "
            ;;
    esac
done 

if [ -n "$MACH" ]; then 
    export MACHINE=$MACH
    if [ $MACHINE == "ls6" ]; then 
        echo "Installing EQquasi on Lonestar6 at TACC ... ..."
        
        echo "Loading netcdf and mumps modules ... ..."
        module load netcdf/4.6.2 mumps
        ml
        
        echo "NETCDF INC and LIB PATH"
        echo $TACC_NETCDF_INC
        echo $TACC_NETCDF_LIB
        
        echo "TACC_MUMPS INC and LIB PATH"
        echo $TACC_MUMPS_INC
        echo $TACC_MUMPS_LIB
        
    elif [ $MACHINE == "ubuntu" ]; then 
        echo "Installing EQquasi on Ubuntu 22.04 ... ..."
        export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
        ln -sf /usr/lib/x86_64-linux-gnu/libblas.so.3 /usr/lib/x86_64-linux-gnu/libblas.so
        ln -sf /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/liblapack.so
    elif [ $MACHINE == "local"]; then 
        MUMPS_LIB_DIR="./mumps/build/local/lib"
        libNames=("libdmumps.a" "libmumps_common.a" "libpord.a" "libsmumps.a")
        all_exist=TRUE
        for file in "${libNames[@]}"; do
            if [! -f "$MUMPS_LIB_DIR/$file" ]; then
                all_exist=FALSE
                break
            fi 
        done 
        
        if $all_exist; then 
            echo "MUMPS have been installed under ./mumps/build/local ..."
        else
            echo "Installing a local copy of MUMPS ..."
            chmod 755 install.mumps.sh
            ./install.mumps.sh
        fi
    fi 
    
    if [ -n "$CONFIG" ]; then 
        echo "Simply configure EQquasi without installation ... ..."
    else
        cd src
        make
        cd ..
        mkdir bin
        mv src/eqquasi bin
    fi

    export EQQUASIROOT=$(pwd)
    export PATH=$(pwd)/bin:$PATH
    export PATH=$(pwd)/scripts:$PATH
    
    echo EQQUASIROOT
    echo PATH 
    
    chmod -R 755 scripts
fi

export EQQUASIROOT=$(pwd)
export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH

echo EQQUASIROOT
echo PATH 

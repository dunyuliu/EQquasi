#! /bin/bash 

# This shell script is to simply set up environments for EQquasi.
# Currently, the machines supported are:
#   ls6:    Lonestar6 at TACC
#   ubuntu: Ubuntu 22.04

# Usage: env_config [-h] [-m] Machine-name

while getopts "hm:" OPTION; do
    case $OPTION in
        m)
            MACH=$OPTARG
            ;;
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

module load netcdf
ml

export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH
export ECCIROOT=$(pwd)
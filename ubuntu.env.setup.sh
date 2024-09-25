#! /bin/bash

# This script contains commands to install 
#   necessary packages for EQquasi on Ubuntu.
apt-get update
apt-get install git vim openmpi-bin make cmake
apt-get install python3 python3-pip  
apt-get install libmumps-dev python3-matplotlib
apt-get install libnetcdf-dev libnetcdff-dev
pip install xarray imageio pdf2image numpy==1.26.4 netCDF4

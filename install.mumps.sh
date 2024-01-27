#! /bin/bash 

# install mumps 
echo 'Installing a local copy of mumps through scivision/mumps ... ...'

rm -rf mumps
git clone https://github.com/scivision/mumps.git
cd mumps
cmake -B build -DMUMPS_UPSTREAM_VERSION=5.6.2
cmake --build build
cd ..

echo 'mumps is installed ...'

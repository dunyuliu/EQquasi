#! /bin/bash 

echo 'Installing a local copy of MUMPS through scivision/mumps ...'
echo 'Deleting old copy ...'
rm -rf mumps

git clone https://github.com/scivision/mumps.git
cd mumps
cmake -B build # use the latest MUMPS.
# cmake -B build -DMUMPS_UPSTREAM_VERSION=5.6.2
cmake --build build
cmake --install build
cd ..

echo 'MUMPS should be installed now ...'
echo 'Libraries are under /mumps/build/local/lib ...'
echo 'Include files are under /mumps/build/local/include ...'

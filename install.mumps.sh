#! /bin/bash 

# install mumps 
echo 'Installing a local copy of mumps ... ...'

rm -rf mumps
git clone https://github.com/scivision/mumps.git
cd mumps
cmake -B build
cmake --build build
cd ..

echo 'mumps is installed ...'

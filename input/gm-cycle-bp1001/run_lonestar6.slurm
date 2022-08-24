#!/bin/bash
#SBATCH -J v1.3.2_bp1001-1000            # job name
#SBATCH -o a.eqquasi.log.o%j        # output and error file name (%j expands to jobID)
#SBATCH -N 1                # number of nodes requested
#SBATCH -n 30               # total number of mpi tasks requested
#SBATCH -p normal      # queue (partition) --small (2), normal(3-512), development, etc.
#SBATCH -t 00:10:00         # run time (hh:mm:ss) - 1.5 hours
#SBATCH -A EAR22013
#SBATCH --mail-user=dliu@ig.utexas.edu
##SBATCH --mail-type=begin   # email me when the job starts
##SBATCH --mail-type=end     # email me when the job finishes

#module load remora
#
for i in $(seq 1 1);
do 
  echo "Cycle #:" $i
  echo "The" $i "cycle started at" >> cycles.log
  date >> cycles.log
  echo $i >> currentcycle.txt
  ibrun eqquasi
  mkdir "./Q$((i-1))"
  mv cplot* fltst* srfst* tdyna* global* p1out* p2out* "./Q$((i-1))"
  cp "./Q$((i-1))/cplot_EQquasi.txt" "./"
  rm currentcycle.txt
  echo "The" $i "cycle ended at" >> cycles.log
done

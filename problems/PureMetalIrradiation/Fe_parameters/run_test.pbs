#!/bin/bash
#######################################
# Specify nodes, processors per node,
# and maximum running time
#######################################
#PBS -N ion

#PBS -l select=1:ncpus=24:mpiprocs=24
#PBS -l walltime=72:00:00
#PBS -q general 
#PBS -P ldrd
##PBS -m e
##PBS -M miaomiaojin1992@gmail.com
###PBS -o 'qsub1.out'
###PBS -e 'qsub1.err'
#######################################
# Enter directory and set PATH
#######################################
#export LD_LIBRARY_PATH=/opt/openmpi-1.6.5/lib:$LD_LIBRARY_PATH
#PATH=$PBS_O_PATH

module load use.moose moose-dev-gcc
#source /home/jinmiao/.bashrc
cd $PBS_O_WORKDIR 
#mpirun /home/jinmiao/inl_projects/falcon/grime-dev/grime-dev-opt -i 042816_golubov.i
#mpirun /home/jinmiao/inl_projects/falcon/grime-dev/grime-dev-opt -i 060316_CD.i
#mpirun /home/jinmiao/inl_projects/falcon/grime-dev/grime-dev-opt -i 042316_jourdan.i
echo "Start: `date`"
mpirun /home/jinmiao/inl_projects/falcon/grime-dev/grime-dev-opt -i test.i
echo "End: `date`"

#mpiexec /home/jinmiao/projects/lammps-16Feb16/src/lmp_mpi < in.cnt_Fe

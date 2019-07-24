#!/bin/bash --login
#PBS -N ozone
#PBS -V
#PBS -l select=1:ncpus=4
#PBS -l walltime=48:00:00
#PBS -A sc007

module load anaconda/python3
export OMP_NUM_THREADS=4

cd $PBS_O_WORKDIR

/lustre/home/sc007/mohammed/ozone/ozone_absorption_for_cirrus.py

#!/bin/bash

make clean
make
qsub job-knapsack
qstat -u omers88

# module swap PrgEnv-pgi PrgEnv-cray
# module load bupc
# qsub -I -V -q interactive -l mppwidth=96
# (wait for node)
# make knapsack && aprun -n 4 -N 1 ./knapsack
# make && aprun -n 4 -N 1 ./knapsack_orig
# make && aprun -n 1 -N 1 ./serial
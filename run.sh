#!/bin/bash

make clean
make
qsub job-knapsack
qstat -u omers88
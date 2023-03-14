#!/bin/bash


cp hubb.andpar_bak hubb.andpar
mpif90 ed_dmft_parallel_frequencies_corrected_and_commented_betterSampling.f -llapack -o run.x  -ffixed-line-length-0  #-fcheck=all
mpirun -np 36 --oversubscribe ./run.x > run.out

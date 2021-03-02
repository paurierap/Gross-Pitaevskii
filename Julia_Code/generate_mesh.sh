#!/bin/bash -l
#SBATCH -J gross-pita
# Define the time allocation you use
#SBATCH -A 2020-72
# 10-seconds wall-clock time will be given to this job
#SBATCH -t 00:00:10
# load Julia compiler.
module add julialang/1.0.1
# Run program
julia Generate_Mesh1D.jl

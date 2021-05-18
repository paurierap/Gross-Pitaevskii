#!/bin/bash -l
#SBATCH -J gross-pita
# Define the time allocation you use
#SBATCH -A 2020-72
# 12-hour  wall-clock time will be given to this job
#SBATCH -t 24:00:00
# Request mail when job starts and ends
#SBATCH --mail-type=END
#SBATCH -C Haswell
# load Julia compiler.
module load julialang/1.5.4-haswell
# Run program
julia main.jl

#!/bin/bash -l
#SBATCH -J gross-pita
# Define the time allocation you use
#SBATCH -A 2020-72
# 12-hour  wall-clock time will be given to this job
#SBATCH -t 24:00:00
# Request mail when job starts and ends
#SBATCH --mail-type=END
# load Julia compiler.
module add julialang/1.0.1
# Run program
julia main.jl

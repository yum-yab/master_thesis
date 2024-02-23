#!/bin/bash
#SBATCH --job-name=generate_ivt_fields
#SBATCH --account=kv0728
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=05:00:00
#SBATCH --partition=compute
#SBATCH --output=output.log
#SBATCH --error=error.log

module load netcdf-c
module load cdo
module load julia/1.7.0-gcc-11.2.0

# it is assumend the env was instantiated correctly
julia --project=$PWD --threads=auto ./script.jl

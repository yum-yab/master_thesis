#!/bin/bash
#SBATCH --job-name=preprocess_ivt_data
#SBATCH --account=kv0728
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=256
#SBATCH --time=05:00:00
#SBATCH --partition=compute
#SBATCH --output=output_xarray_loading.log
#SBATCH --error=error_xarray_loading.log
#SBATCH --mem=80G

module load netcdf-c
module load cdo
module load julia/1.7.0-gcc-11.2.0
module load python3

source env/bin/activate

julia --project --threads=auto script.jl


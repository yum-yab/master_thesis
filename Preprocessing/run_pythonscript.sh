#!/bin/bash
#SBATCH --job-name=test_python_xarray_pure
#SBATCH --account=kv0728
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=256
#SBATCH --time=01:00:00
#SBATCH --partition=compute
#SBATCH --output=output_xr_impl.log
#SBATCH --error=error_xr_impl.log
#SBATCH --mem=80G
module load netcdf-c
module load python3

source env/bin/activate

python3 vertical_integration_xarray.py 

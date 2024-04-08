#!/bin/bash
#SBATCH --job-name=test_python_xarray_pure
#SBATCH --account=kv0728
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=256
#SBATCH --time=01:00:00
#SBATCH --partition=compute
#SBATCH --output=output_full_member_run.log
#SBATCH --error=error_full_member_run.log
#SBATCH --mem=100G
module load netcdf-c
module load python3

source env/bin/activate

python3 vertical_integration_xarray.py $1 $2

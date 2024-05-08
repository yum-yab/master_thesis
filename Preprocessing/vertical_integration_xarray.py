import xarray as xr 
import os
import sys
import glob
import numpy as np
from dask.distributed import Client
import time
import multiprocessing

# In order to run, this script needs the following 3rd party libraries
#
# Requirements:
#   xarray[complete]
#   dask[complete]
#   numpy
#   multiprocessing



def integrate_trapz(array, x):
    return np.trapz(array, x, axis=-1) / 9.806

def integrate_over_data_array(y_values, x_values, integration_dim = 'lev'):
    return xr.apply_ufunc(
        integrate_trapz,
        y_values,
        x_values,
        input_core_dims=[[integration_dim], [integration_dim]],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[float]  # Specify the output data type
     )



def generate_ivt_calculation(source_paths, target_path, chunks=dict(time=256, lev= 47, lat= 96, lon= 192)):

    start_time = time.perf_counter()

    ds = xr.open_mfdataset(source_paths, compat="override", decode_cf=False, parallel=True, chunks=chunks)
    # ds = xr.open_mfdataset(paths, compat="override",  decode_cf=False, chunks=chunks)

    opened_datasets = time.perf_counter()
    print(f"Time used for loading dataset: {opened_datasets - start_time} s")

    ds.coords["lon"] = (ds.coords["lon"] + 180) % 360 - 180
        
    relevant_subset = ds.sortby(ds.lon).sel(lon=slice(-90, 40), lat=slice(20, 80))

    relevant_subset["plevs"] = relevant_subset["ap"] + relevant_subset["b"] * relevant_subset["ps"]

    eastward_product = relevant_subset["hus"] * relevant_subset["ua"]

    northward_product = relevant_subset["hus"] * relevant_subset["va"]

    eastward_integral = integrate_over_data_array(eastward_product, relevant_subset.plevs)

    northward_integral = integrate_over_data_array(northward_product, relevant_subset.plevs)

    ivt_norm = np.sqrt(eastward_integral**2 + northward_integral**2)

    new_ds = xr.Dataset(coords={coord: relevant_subset[coord] for coord in relevant_subset.coords if coord != 'lev'})

    new_ds["time"] = relevant_subset.time

    new_ds.attrs = relevant_subset.attrs

    new_ds.attrs["variable_id"] = "ivt"

    new_ds["ivt_east_component"] = eastward_integral
    new_ds["ivt_east_component"].attrs.update(
        long_name="Integrated Water Vapor Transport Eastward Component",
        units="kg m^-1 s^-1",
        comments="Eastward water vapor flux. Calculated by integrating over hus * ua fields. This field is not normalized.",
        history="Genrated from MPI-M CMIP6 simulation by integrating with xarray+dask. See https://github.com/yum-yab/master_thesis"
    )
    new_ds["ivt_north_component"] = northward_integral
    new_ds["ivt_north_component"].attrs.update(
        long_name="Integrated Water Vapor Transport Northward Component",
        units= "kg m^-1 s^-1",
        comments= "Northward water vapor flux. Calculated by integrating over hus * va fields. This field is not normalized.",
        history="Genrated from MPI-M CMIP6 simulation by integrating with xarray+dask. See https://github.com/yum-yab/master_thesis"
    )
    new_ds["ivt"] = ivt_norm
    new_ds["ivt"].attrs.update(
        long_name="Integrated Water Vapor Transport normalized",
        units="kg m^-1 s^-1",
        comment="Euclidian norm of the water vapor flux vector. Calculated by (ivteast^2+ivtnorth^2)^0.5. This field is not normalized.",
        history="Genrated from MPI-M CMIP6 simulation by integrating with xarray+dask. See https://github.com/yum-yab/master_thesis"
    )

    all_lazy_calcs = time.perf_counter()
    print(f"Time used for lazy calculations: {all_lazy_calcs - opened_datasets}")

    new_ds.to_netcdf(target_path)

    very_end = time.perf_counter()

    print(f"Time used for actual calculation: {very_end - all_lazy_calcs} s")
    print(f"Time used all in all: {(very_end - start_time)/60} min")


def get_latest_version(base_path):

    all_versions = [de.name for de in os.scandir(base_path) if de.is_dir()]
    all_versions.sort()

    return all_versions[-1]

def timestamp_from_nc_file(nc_file_name):
    """
    Extracts and returns a timestamp string from a netCDF file name.
    """
    return nc_file_name.split("_")[-1].split(".")[0]

def get_files_per_timescope(member_base, field_ids, time_res_id = "6hrLev"):

    sorted_files = []

    for field_id in field_ids:
        
        path = os.path.join(member_base, time_res_id, field_id, "gn")
        
        latest_version = get_latest_version(path)
        
        version_path = os.path.join(path, latest_version)

        nc_files = glob.glob(os.path.join(version_path, "*.nc"))

        nc_files.sort()

        sorted_files.append(nc_files)

    for files in zip(*sorted_files):

        files_for_timescope = list(files)

        yield (timestamp_from_nc_file(files_for_timescope[0]), files_for_timescope)


    


    
if __name__ == "__main__":
    
    ncpu = multiprocessing.cpu_count()
    processes = True
    nworker = 16
    threads = ncpu // nworker
    mem_limit_per_worker="24GB"
    print(
        f"Number of CPUs: {ncpu}, number of threads: {threads}, number of workers: {nworker}, processes: {processes}, mem_limit_per_worker: {mem_limit_per_worker}",
    )
    
    field_ids = ["hus", "ua", "va"]

    arguments = sys.argv[1:]

    member_base_path = arguments[0]

    scenario_path, member_id = os.path.split(member_base_path)

    _, scenario_name = os.path.split(scenario_path)

    target_base = arguments[1]

    if len(arguments) == 3:
        time_res_id = arguments[2]
    else:
        time_res_id = "6hrLev"
    
    with Client(processes=processes, threads_per_worker=threads, n_workers=nworker, memory_limit=mem_limit_per_worker) as client:
        # dask.config.config.get('distributed').get('dashboard').update({'link':'{JUPYTERHUB_SERVICE_PREFIX}/proxy/{port}/status'})

        for timestamp, field_files in get_files_per_timescope(member_base_path, field_ids, time_res_id = time_res_id):

            target_path = os.path.join(target_base, scenario_name, member_id)

            os.makedirs(target_path, exist_ok=True)

            target_file = os.path.join(target_path, f"ivt_{scenario_name}_{member_id}_{timestamp}.nc")
            
            print(f"Generates ivt for files {field_files} and target {target_file}")
            generate_ivt_calculation(field_files, target_file, chunks=dict(time=128, lev= 47, lat= 96, lon= 192))

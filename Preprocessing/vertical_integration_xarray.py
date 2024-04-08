import xarray as xr 
import numpy as np
from dask.distributed import LocalCluster
import time
import multiprocessing





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

def generate_ivt_calculation(source_paths, target_path, chunks=dict(time=64, lev= 47, lat= 96, lon= 192)):

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

    new_ds.attrs = relevant_subset.attrs

    new_ds["eastward_integral"] = eastward_integral
    new_ds["northward_integral"] = northward_integral
    new_ds["ivt"] = ivt_norm

    all_lazy_calcs = time.perf_counter()
    print(f"Time used for lazy calculations: {all_lazy_calcs - opened_datasets}")

    new_ds.to_netcdf(target_path)

    very_end = time.perf_counter()

    print(f"Time used for actual calculation: {very_end - all_lazy_calcs} s")
    print(f"Time used all in all: {(very_end - start_time)/60} min")

    
    
if __name__ == "__main__":
    
    ncpu = multiprocessing.cpu_count()
    processes = True
    nworker = 16
    threads = ncpu // nworker
    mem_limit_per_worker="16GB"
    print(
        f"Number of CPUs: {ncpu}, number of threads: {threads}, number of workers: {nworker}, processes: {processes}, mem_limit_per_worker: {mem_limit_per_worker}",
    )

    
    with Client(processes=processes, threads_per_worker=threads, n_workers=nworker, memory_limit=mem_limit_per_worker) as client:
        dask.config.config.get('distributed').get('dashboard').update({'link':'{JUPYTERHUB_SERVICE_PREFIX}/proxy/{port}/status'})
        paths = [
            "/work/ik1017/CMIP6/data/CMIP6/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp585/r10i1p1f1/6hrLev/hus/gn/v20190710/hus_6hrLev_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501010600-203501010000.nc", 
            "/work/ik1017/CMIP6/data/CMIP6/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp585/r10i1p1f1/6hrLev/va/gn/v20190815/va_6hrLev_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501010600-203501010000.nc", 
            "/work/ik1017/CMIP6/data/CMIP6/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp585/r10i1p1f1/6hrLev/ua/gn/v20190815/ua_6hrLev_MPI-ESM1-2-LR_ssp585_r10i1p1f1_gn_201501010600-203501010000.nc"
        ]


        target_path = "/scratch/b/b382641/test_script_ssp585_r10i1p1f1.nc.nc"

        generate_ivt_calculation(paths, target_path)
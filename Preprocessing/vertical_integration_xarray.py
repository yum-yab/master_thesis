import xarray as xr 
import numpy as np


paths = ["/home/denis/Documents/Uni/Master/MA/preprocessing/sample_data/sample_hus_dataset_200_timesteps.nc", "/home/denis/Documents/Uni/Master/MA/preprocessing/sample_data/sample_va_dataset_200_timesteps.nc", "/home/denis/Documents/Uni/Master/MA/preprocessing/sample_data/sample_ua_dataset_200_timesteps.nc"]


target_path = "./test_xarray_impl.nc"
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
    

ds = xr.open_mfdataset(paths, compat="override", parallel=True, chunks=dict(time=256, lev= 47, lat= 96, lon= 192))

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

new_ds.to_netcdf(target_path)

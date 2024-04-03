import xarray as xr 


paths = ["/home/denis/Documents/Uni/Master/MA/Preprocessing/sample_data/sample_hus_dataset_200_timesteps.nc", "/home/denis/Documents/Uni/Master/MA/Preprocessing/sample_data/sample_va_dataset_200_timesteps.nc", "/home/denis/Documents/Uni/Master/MA/Preprocessing/sample_data/sample_ua_dataset_200_timesteps.nc"]

ds = xr.open_mfdataset(paths, compat="override", parallel=True, chunks=dict(time=256, lev= 47, lat= 96, lon= 192))

ds.coords["lon"] = (ds.coords["lon"] + 180) % 360 - 180
    

relevant_subset = ds.sortby(ds.lon).sel(lon=slice(-90, 40), lat=slice(20, 80))

relevant_subset["plevs"] = ds["ap"] + ds["b"] * ds["ps"]

eastward_product = ds["hus"] * ds["ua"]

northward_product = ds["hus"] * ds["va"]

eastward_integral = eastward_product.integrate("")

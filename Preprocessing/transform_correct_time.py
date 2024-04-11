import os
import xarray as xr
import pandas as pd

def transform_and_update_metadata(ds):
    # Transform the time coordinate
    reference_date = pd.Timestamp('1850-01-01')
    ds['time'] = pd.to_datetime(ds['time'].values, unit='D', origin=reference_date)
    # Update metadata for time, removing conflicting keys if they exist
    time_metadata = {
        "axis": "T",
        "long_name": "time",
        "standard_name": "time"
    }
    for key in ['units', 'calendar']:
        ds['time'].attrs.pop(key, None)
    ds['time'].attrs.update(time_metadata)
    ds['time'].encoding['dtype'] = 'float64'
    ds['time'].encoding['units'] = "days since 1850-1-1 00:00:00"
    ds['time'].encoding['calendar'] = "proleptic_gregorian"

    # Update metadata for lon and lat, removing conflicting keys if they exist
    lon_metadata = {
        "axis": "X",
        "long_name": "Longitude",
        "standard_name": "longitude"
    }
    lat_metadata = {
        "axis": "Y",
        "long_name": "Latitude",
        "standard_name": "latitude"
    }
    for key in ['units']:
        ds['lon'].attrs.pop(key, None)
        ds['lat'].attrs.pop(key, None)
    ds['lon'].attrs.update(lon_metadata)
    ds['lat'].attrs.update(lat_metadata)
    ds['lon'].encoding['units'] = "degrees_east"
    ds['lat'].encoding['units'] = "degrees_north"

    return ds

def process_files(rootdir, rootdir_new):
    for dirpath, dirnames, filenames in os.walk(rootdir):
        for filename in filenames:
            if filename.endswith('.nc'):

                # Construct the old and new file paths
                old_file_path = os.path.join(dirpath, filename)
                new_file_path = old_file_path.replace(rootdir, rootdir_new)

                print(f"Handling file {old_file_path}")
                # Ensure the new directory exists
                os.makedirs(os.path.dirname(new_file_path), exist_ok=True)

                # Load, transform, and save the dataset
                ds = xr.open_dataset(old_file_path)
                ds = transform_and_update_metadata(ds)
                ds.to_netcdf(new_file_path)
                ds.close()

# Example usage
rootdir = '/home/denis/workspace/data/ivt_fields_v1'
rootdir_new = '/home/denis/workspace/data/ivt_fields_correct_time'
process_files(rootdir, rootdir_new)

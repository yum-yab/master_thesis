import os
import xarray as xr
import pandas as pd

def transform_and_update_metadata(ds):
    """Transform the time coordinate and update metadata for time, lon, and lat."""
    # Transform the time coordinate
    reference_date = pd.Timestamp('1850-01-01')
    ds['time'] = pd.to_datetime(ds['time'].values, unit='D', origin=reference_date)

    # Metadata for time, longitude, and latitude
    time_metadata = {
        "units": "days since 1850-1-1 00:00:00",
        "calendar": "proleptic_gregorian",
        "axis": "T",
        "long_name": "time",
        "standard_name": "time"
    }
    lon_metadata = {
        "units": "degrees_east",
        "axis": "X",
        "long_name": "Longitude",
        "standard_name": "longitude"
    }
    lat_metadata = {
        "units": "degrees_north",
        "axis": "Y",
        "long_name": "Latitude",
        "standard_name": "latitude"
    }

    # Update metadata
    ds['time'].attrs.update(time_metadata)
    ds['lon'].attrs.update(lon_metadata)
    ds['lat'].attrs.update(lat_metadata)

    return ds

def process_files(rootdir, rootdir_new):
    for dirpath, dirnames, filenames in os.walk(rootdir):
        for filename in filenames:
            if filename.endswith('.nc'):
                # Construct the old and new file paths
                old_file_path = os.path.join(dirpath, filename)
                new_file_path = old_file_path.replace(rootdir, rootdir_new)

                # Ensure the new directory exists
                os.makedirs(os.path.dirname(new_file_path), exist_ok=True)

                # Load, transform, and save the dataset
                ds = xr.open_dataset(old_file_path)
                ds = transform_and_update_metadata(ds)
                ds.to_netcdf(new_file_path)
                ds.close()

# Example usage
rootdir = '/path/to/rootdir'
rootdir_new = '/path/to/rootdir_new'
process_files(rootdir, rootdir_new)

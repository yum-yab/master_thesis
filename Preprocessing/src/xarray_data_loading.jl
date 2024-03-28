
module XarrayDataLoading
  
using Distributed

using PythonCall
export parallel_loading_of_datasets


function parallel_loading_of_datasets(id_to_file_mapping::Dict{String, String})::Dict{String, Array{<: AbstractFloat}}
  
  arguments = zip([(id, path) for (id, path) in id_to_file_mapping]...) 
  
  results = pmap(arguments...) do id, path
    data = load_one_dataset_xarray(path, id)
    return id => data
  end
  
  return Dict(results)
end
function load_one_dataset_xarray(path::String, id::String)::AbstractArray{<: AbstractFloat}
    
    println("Started loading data $id from $path")
    xr = pyimport("xarray")
    ds = xr.open_dataset(path, chunks=pydict(lon = 192, lat = 96, lev = 47, time= 256))
    ds.coords["lon"] = (ds.coords["lon"] + 180) % 360 - 180
    relevant_subset = ds.sortby(ds.lon).sel(lon=pyslice(-90, 40), lat=pyslice(20, 80))
    
    data = relevant_subset[id].to_numpy() 
  
    ds.close()
  
    return pyconvert(Array, data)
    
  end
end



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
    
    println("Hello from worker $id")
    xr = pyimport("xarray")
    ds = xr.open_dataset(path, chunks=pydict(lon = 192, lat = 96, lev = 47, time= 256))
  
    ds.assign_coords(Dict("lon" => (ds.lon.values .+ 180) .% 360 .- 180))
    relevant_subset = ds.sortby(ds.lon).sel(lon=pyslice(-90, 40), lat=pyslice(20, 80))
    
    println("Actually loading the data")
    @time data = relevant_subset[id].to_numpy() 
  
    ds.close()
  
    return pyconvert(Array, data)
    
  end
end


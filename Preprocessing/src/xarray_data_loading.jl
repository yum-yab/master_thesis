
module XarrayDataLoading
  
using Distributed
using PythonCall
export parallel_loading_of_datasets

function shared_loading_of_datasets(id_to_file_mapping::Dict{String, String}, dims)::Dict{String, Array{<: AbstractFloat}}
  
  hus_array = SharedArray{Float32}(dims...)

  results = pmap(arguments...) do id, path
    data = load_one_dataset_xarray(path, id)
    return id => data
  end
  
  return Dict(results)
end

function parallel_loading_of_datasets(id_to_file_mapping::Dict{String, String})::Dict{String, Array{<: AbstractFloat}}
  
  arguments = zip([(id, path) for (id, path) in id_to_file_mapping]...) 
  
  results = pmap(arguments...) do id, path
    data = load_one_dataset_xarray(path, id)
    return id => data
  end

  for pair in results
    println("Size of compressed $(pair[1]) data: $(sizeof(pair[2]) / 10^6) MB")
  end
  
  return Dict([p[1] => decompress_data(p[2]) for p in results])
end

function iterative_loading_of_datasets(id_to_file_mapping::Dict{String, String})

  return Dict([id => load_one_dataset_xarray(path, id) for (id, path) in id_to_file_mapping])
  
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

function multithreaded_loading_of_datasets(id_to_file_mapping::Dict{String, String})::Dict{String, Array{<: AbstractFloat}}
  
  res = Dict{String, AbstractArray}()

  arguments = [(id, path) for (id, path) in id_to_file_mapping] 
  Threads.@threads for i in eachindex(arguments)
    id = arguments[i][1]
    path = arguments[i][2]
    res[id] = load_one_dataset_xarray(path, id)
  end  
  return res 
end
end


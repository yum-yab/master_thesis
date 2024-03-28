
module XarrayDataLoading
  

using Distributed
export remap_and_load_data_with_xarray
@everywhere using PythonCall

function remap_and_load_data_with_xarray(files::Vector{String})::Dict{String, Array{<: AbstractFloat}}
  
  xr = pyimport("xarray")
  println("Opening the dataset")
  @time ds = xr.open_mfdataset(files, chunks=py"dict(time=256, lev= 47, lat= 96, lon= 192)", compat="override")
  
  ds.assign_coords(Dict("lon" => (ds.lon.values .+ 180) .% 360 .- 180))
  relevant_subset = ds.sortby(ds.lon).sel(lon=pyslice(-90, 40), lat=pyslice(20, 80))
  
  println("Actually loading the data")
  @time result_dict = Dict([id => relevant_subset[id].to_numpy() for id in ["ap", "b", "ps", "hus", "ua", "va"]])

  ds.close()

  return result_dict
end

function parallel_loading_of_datasets(id_to_file_mapping::Dict{String, String})::Dict{String, Array{<: AbstractFloat}}
  result_dict = Dict{String, AbstractArray}()
  arguments = zip([(id, path) for (id, path) in id_to_file_mapping]...) 
  addprocs(length(id_to_file_mapping))
  barrier = ReentrantLock()
  pmap(arguments...) do id, path
    data = load_one_dataset_xarray(path, id)
    lock(barrier)
    result_dict[id] = data
    unlock(barrier)
  end
  return result_dict
end

@everywhere function load_one_dataset_xarray(path::String, id::String)::AbstractArray{<: AbstractFloat}
  
  println("Hello from worker $id")
  @everywhere xr = pyimport("xarray")
  ds = xr.open_dataset(path, chunks=Dict("lon" => 192, "lat" => 96, "lev" => 47, "time" => 256))

  ds.assign_coords(Dict("lon" => (ds.lon.values .+ 180) .% 360 .- 180))
  relevant_subset = ds.sortby(ds.lon).sel(lon=pyslice(-90, 40), lat=pyslice(20, 80))
  
  println("Actually loading the data")
  @time data = relevant_subset[id].to_numpy() 

  ds.close()

  return data
  
end
end


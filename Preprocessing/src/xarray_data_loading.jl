
module XarrayDataLoading
  
ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")

using PyCall 
xr = PyNULL()

export remap_and_load_data_with_xarray

function __init__()
  copy!(xr, pyimport("xarray"))
end

function remap_and_load_data_with_xarray(files::Vector{String})::Dict{String, Array{<: AbstractFloat}}
  
  println("Opening the dataset")
  ds = xr.open_mfdataset(files, chunks=py"dict(time=256, lev= 47, lat= 96, lon= 192)")
  
  println("Remapping longitude coordinates")
  ds.coords["lon"] = (ds.coords["lon"] + 180) % 360 - 180
  relevant_subset = ds.sortby(ds.lon).sel(lon=py"slice(-90, 40)", lat=py"slice(20, 80)")
  
  println("Actually loading the data")
  result_dict = Dict([id => relevant_subset[id].to_numpy() for id in ["ap", "b", "ps", "hus", "ua", "va"]])

  ds.close()

  return result_dict
end
end


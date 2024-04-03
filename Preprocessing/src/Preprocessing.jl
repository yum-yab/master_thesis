
module Preprocessing 

using NCDatasets
using DataStructures
using NetCDF

include("IVT.jl")
using .IVT

include("data_loading.jl")
using .DataLoading

include("xarray_data_loading.jl")
using .XarrayDataLoading

export generate_ivt_field, write_ivt_dataset, IVT, DataLoading, generate_ivt_field_xarray_loading, XarrayDataLoading



function generate_ivt_field(id_to_file_mapping::Dict{String, String}, T::Type = Float64)


  println("Time used for loading the data: ")
  @time begin
    hus_data = ncread(id_to_file_mapping["hus"], "hus")
    ua_data = ncread(id_to_file_mapping["ua"], "ua")
    va_data = ncread(id_to_file_mapping["va"], "va")
    ps_data = ncread(id_to_file_mapping["ps"], "ps")
    
    (lon_size, lat_size, _, time_size) = size(hus_data)

    # these variables are used for calculation of pressure levels at each specific lat, lon, time coordinate: p = ap + b * ps
    # using netcdf lib here since it will be quite fast
    ap = ncread(id_to_file_mapping["hus"], "ap")
    b = ncread(id_to_file_mapping["hus"], "b")
  end

  result_data_eastwards = zeros(T, lon_size, lat_size, time_size)
  result_data_northwards = zeros(T, lon_size, lat_size, time_size)
  result_data_norm = zeros(T, lon_size, lat_size, time_size)
  println("Time used for calculating the IVT field: ")
  @time begin
    Threads.@threads for time in 1:time_size
      for lat in 1:lat_size
        for lon in 1:lon_size

          ps = ps_data[lon, lat, time]
          pressure_levels = ap + b * ps

          result = IVT.ivt_of_column(
            pressure_levels, 
            hus_data[lon, lat, :, time], 
            ua_data[lon, lat, :, time], 
            va_data[lon, lat, :, time]
          )

          result_data_eastwards[lon, lat, time] = result.eastward_integral
          result_data_northwards[lon, lat, time] = result.northward_integral
          result_data_norm[lon, lat, time] = sqrt(result.northward_integral^2 + result.eastward_integral^2)



        end
      end
    end

  end

  return result_data_eastwards, result_data_northwards, result_data_norm 
end


function generate_ivt_field(id_to_data_mapping::Dict{String,<: Array}, T::Type = Float64)
  return generate_ivt_field(id_to_data_mapping["hus"], id_to_data_mapping["ua"], id_to_data_mapping["va"], id_to_data_mapping["ps"], id_to_data_mapping["ap"], id_to_data_mapping["b"], T) 
end

function generate_ivt_field(hus_data::Array{<: AbstractFloat, 4}, ua_data::Array{<: AbstractFloat, 4}, va_data::Array{<: AbstractFloat, 4}, ps_data::Array{<: AbstractFloat, 3}, ap_data::Array{<: AbstractFloat, 1}, b_data::Array{<: AbstractFloat, 1}, T::Type = Float64)
  
  (time_size, _, lat_size, lon_size) = size(hus_data)
  result_data_eastwards = zeros(T, lon_size, lat_size, time_size)
  result_data_northwards = zeros(T, lon_size, lat_size, time_size)
  result_data_norm = zeros(T, lon_size, lat_size, time_size)
  println("Time used for calculating the IVT field: ")
  @time begin
    Threads.@threads for time in 1:time_size
      for lat in 1:lat_size
        for lon in 1:lon_size

          ps = ps_data[lon, lat, time]
          pressure_levels = ap_data + b_data * ps

          result = IVT.ivt_of_column(
            pressure_levels, 
            hus_data[time, :, lat, lon], 
            ua_data[time, :, lat, lon], 
            va_data[time, :, lat, lon]
          )

          result_data_eastwards[lon, lat, time] = result.eastward_integral
          result_data_northwards[lon, lat, time] = result.northward_integral
          result_data_norm[lon, lat, time] = sqrt(result.northward_integral^2 + result.eastward_integral^2)



        end
      end
    end

  end

  return result_data_eastwards, result_data_northwards, result_data_norm 
end

function generate_ivt_field(id_to_file_mapping::Dict{String, String}, geo_bnds::GeographicBounds, T::Type = Float64)


  println("Time used for loading the data: ")
  @time begin
    hus_data = load_variable_data_in_bounds(Float32, id_to_file_mapping["hus"], "hus", geo_bnds, :, :)
    ua_data = load_variable_data_in_bounds(Float32, id_to_file_mapping["ua"], "ua", geo_bnds, :, :)
    va_data = load_variable_data_in_bounds(Float32, id_to_file_mapping["va"], "va", geo_bnds, :, :)
    ps_data = load_variable_data_in_bounds(Float32, id_to_file_mapping["ps"], "ps", geo_bnds, :)

    lon_size = size(hus_data, 1)
    lat_size = size(hus_data, 2)

    # these variables are used for calculation of pressure levels at each specific lat, lon, time coordinate: p = ap + b * ps
    # using netcdf lib here since it will be quite fast
    ap = ncread(id_to_file_mapping["hus"], "ap")
    b = ncread(id_to_file_mapping["hus"], "b")
    time_size = size(hus_data, 4)

  end

  result_data_eastwards = zeros(T, lon_size, lat_size, time_size)
  result_data_northwards = zeros(T, lon_size, lat_size, time_size)
  result_data_norm = zeros(T, lon_size, lat_size, time_size)
  println("Time used for calculating the IVT field: ")
  @time begin
    Threads.@threads for time in 1:time_size
      for lat in 1:lat_size
        for lon in 1:lon_size

          ps = ps_data[lon, lat, time]
          pressure_levels = ap + b * ps

          result = IVT.ivt_of_column(
            pressure_levels, 
            hus_data[lon, lat, :, time], 
            ua_data[lon, lat, :, time], 
            va_data[lon, lat, :, time]
          )

          result_data_eastwards[lon, lat, time] = result.eastward_integral
          result_data_northwards[lon, lat, time] = result.northward_integral
          result_data_norm[lon, lat, time] = sqrt(result.northward_integral^2 + result.eastward_integral^2)



        end
      end
    end

  end

  return result_data_eastwards, result_data_northwards, result_data_norm 
end

function write_ivt_dataset(old_dataset_attributes::Dict, geo_bnds::GeographicBounds, time_data::Vector, target_path::String, data_eastwards, data_northwards, data_norm)::Nothing
  # copy the data to a new Ordered Dict since a direct pass would only pass the reference to the old atrrib dict -> ERROR!
  new_attrib = OrderedDict()

  for (k, v) in old_dataset_attributes
    new_attrib[k] = v
  end


  # then rewrite some fields
  new_attrib["history"] = "Used an IVT julia script found here: https://github.com/yum-yab/master_thesis/tree/main/preprocessing" * new_attrib["history"]
  new_attrib["title"] = "IVT field of europe and the north-east Atlantic"

  NCDataset(target_path, "c", attrib=new_attrib) do ds
    
    write_dataset_vars(ds, geo_bnds.remaining_lon_values, geo_bnds.remaining_lat_values, time_data, data_eastwards, data_northwards, data_norm)


      return
  end
  return
end
function write_dataset_vars(ds, lon_vals::Vector{<:AbstractFloat}, lat_vals::Vector{<:AbstractFloat}, time_vals::Vector, data_eastwards::Array, data_northwards::Array, data_norm::Array)
  defVar(ds, "lon", lon_vals, ("lon",), attrib=OrderedDict(
    "units" => "degrees east",
    "standard_name" => "longitude"
  ))
  defVar(ds, "lat", lat_vals, ("lat",), attrib=OrderedDict(
    "units" => "degrees north",
    "standard_name" => "latitude"
  ))
  defVar(ds, "time", time_vals, ("time",), attrib=OrderedDict(
    "units" => "days since 1850-1-1 00:00:00",
    "standard_name" => "time"
  ))
  defVar(ds, "ivtnorm", data_norm, ("lon", "lat", "time"), attrib=OrderedDict(
    "units" => "kg/ms",
    "comments" => "Euclidian norm of the water vapor flux vector. Calculated by (ivteast^2+ivtnorth^2)^0.5. This field is not normalized."
  ))
  defVar(ds, "ivteast", data_eastwards, ("lon", "lat", "time"), attrib=OrderedDict(
    "units" => "kg/ms",
    "comments" => "Eastward water vapor flux. Calculated by integrating over hus * ua fields. This field is not normalized."
  ))
  defVar(ds, "ivtnorth", data_northwards, ("lon", "lat", "time"), attrib=OrderedDict(
    "units" => "kg/ms",
    "comments" => "Northward water vapor flux. Calculated by integrating over hus * va fields. This field is not normalized."
  ))
  return
end
function write_ivt_dataset(old_dataset_path::String, target_path::String, data_eastwards, data_northwards, data_norm)::Nothing
  NCDataset(old_dataset_path) do old_ds
    # copy the data to a new Ordered Dict since a direct pass would only pass the reference to the old atrrib dict -> ERROR!
    new_attrib = OrderedDict()

    for (k, v) in old_ds.attrib
      new_attrib[k] = v
    end


    # then rewrite some fields
    new_attrib["history"] = "Used an IVT julia script found here: https://github.com/yum-yab/master_thesis/tree/main/Preprocessing" * new_attrib["history"]
    new_attrib["title"] = "IVT field of europe and the north-east Atlantic"

    NCDataset(target_path, "c", attrib=new_attrib) do ds

      defVar(ds, "lon", old_ds["lon"][:], ("lon",), attrib=OrderedDict(
        "units" => "degrees east",
        "standard_name" => "longitude"
      ))
      defVar(ds, "lat", old_ds["lat"][:], ("lat",), attrib=OrderedDict(
        "units" => "degrees north",
        "standard_name" => "latitude"
      ))
      defVar(ds, "time", old_ds[:time][:], ("time",), attrib=OrderedDict(
        "units" => "days since 1850-1-1 00:00:00",
        "standard_name" => "time"
      ))
      defVar(ds, "ivtnorm", data_norm, ("lon", "lat", "time"), attrib=OrderedDict(
        "units" => "kg/ms",
        "comments" => "Euclidian norm of the water vapor flux vector. Calculated by (ivteast^2+ivtnorth^2)^0.5. This field is not normalized."
      ))
      defVar(ds, "ivteast", data_eastwards, ("lon", "lat", "time"), attrib=OrderedDict(
        "units" => "kg/ms",
        "comments" => "Eastward water vapor flux. Calculated by integrating over hus * ua fields. This field is not normalized."
      ))
      defVar(ds, "ivtnorth", data_northwards, ("lon", "lat", "time"), attrib=OrderedDict(
        "units" => "kg/ms",
        "comments" => "Northward water vapor flux. Calculated by integrating over hus * va fields. This field is not normalized."
      ))


        return
    end
  end
  return
end
end


module preprocessing

  using NCDatasets
  using Distributed
  using NCDatasets
  using DataStructures
  
  include("IVT.jl")
  using .IVT


  export generate_ivt_field
  
  struct GeographicBounds
    """This struct sets the geographic boundaries of the datasets. """
    lon_bounds::Tuple{<:Real, <:Real}
    lon_indices::Tuple{<:Int, <:Int}
    remaining_lon_values::Vector{<:Real}
    lat_bounds::Tuple{<:Real, <:Real}
    lat_indices::Tuple{<:Int, <:Int}
    remaining_lat_values::Vector{<:Real}
  end

  function GeographicBounds(lon_bounds::Tuple{<:Real, <:Real}, lat_bounds::Tuple{<:Real, <:Real}, lon_values::Vector{<:Real}, lat_values::Vector{<:Real})
  
    function get_overflowing_result(value_bounds, values)::Tuple{Tuple{<:Int, <:Int}, Vector{<:Real}}
      
      first_id = findfirst(x -> x >= value_bounds[1], values)
      second_id = findlast(x -> x <= value_bounds[2], values)
      index_bounds = (first_id, second_id)  

      if value_bounds[1] < value_bounds[2]
        remaining_values = values[index_bounds[1]:index_bounds[2]]
      else
        remaining_values = vcat(values[index_bounds[1]:end], values[1:index_bounds[2]]) 
      end
      return index_bounds, remaining_values
    end
    
    (lon_indices, remaining_lon_vals) = get_overflowing_result(lon_bounds, lon_values)
    (lat_indices, remaining_lat_vals) = get_overflowing_result(lat_bounds, lat_values)
    
    return GeographicBounds(lon_bounds, lon_indices, remaining_lon_vals, lat_bounds, lat_indices, remaining_lat_vals)
  end

  function GeographicBounds(lon_bounds::Tuple{<:Real, <:Real}, lat_bounds::Tuple{<:Real, <:Real}, dataset, lon_id = :lon, lat_id = :lat)
    lon_vals = dataset[lon_id][:]
    lat_vals = dataset[lat_id][:]

    return GeographicBounds(lon_bounds, lat_bounds, lon_vals, lat_vals)
  end
  

  
  
"""This function loads data in given geographic bounds. It supports loading values going 'over the end' like the lon rage from 270-40 NOTE: It is expected that the longitude is given in values from 0-360 deg and lat in range from -90:90""" 
function load_data_in_geo_bounds(dataset, field_id::Union{String, Symbol, Missing}, geo_bnds::GeographicBounds, indices...)::Array

  lon_normal_range = geo_bnds.lon_bounds[1] < geo_bnds.lon_bounds[2]
  lat_normal_range = geo_bnds.lat_bounds[1] < geo_bnds.lat_bounds[2]


  if lon_normal_range & lat_normal_range
    return dataset[field_id][geo_bnds.lon_indices[1]:geo_bnds.lon_indices[2], geo_bnds.lat_indices[1]:geo_bnds.lat_indices[2], indices...]
  elseif !lon_normal_range & lat_normal_range
    
    lon_first = dataset[field_id][geo_bnds.lon_indices[1]:end, geo_bnds.lat_indices[1]:geo_bnds.lat_indices[2],indices...]
    lon_second = dataset[field_id][1:geo_bnds.lon_indices[2], geo_bnds.lat_indices[1]:geo_bnds.lat_indices[2],indices...]
    
    return vcat(lon_first, lon_second)
  elseif lon_normal_range & !lat_normal_range
    
    lat_first = dataset[field_id][geo_bnds.lon_indices[1]:geo_bnds.lon_indices[2], geo_bnds.lat_indices[1]:end,indices...]
    lat_second = dataset[field_id][geo_bnds.lon_indices[1]:geo_bnds.lon_indices[2], 1:geo_bnds.lat_indices[2],indices...]
    
    return hcat(lat_first, lat_second)
  else
    # last case is both are over 
    lon_f_lat_f = dataset[field_id][geo_bnds.lon_indices[1]:end, geo_bnds.lat_indices[1]:end,indices...]
    lon_f_lat_s = dataset[field_id][geo_bnds.lon_indices[1]:end, 1:geo_bnds.lat_indices[2],indices...]

    lon_s_lat_f = dataset[field_id][1:geo_bnds.lon_indices[2], geo_bnds.lat_indices[1]:end,indices...]
    lon_s_lat_s = dataset[field_id][1:geo_bnds.lon_indices[2], 1:geo_bnds.lat_indices[2],indices...]

    return vcat(hcat(lon_f_lat_f, lon_f_lat_s), hcat(lon_s_lat_f, lon_s_lat_s))
  end
end
  

  function generate_ivt_field(relevant_fields::Union{Vector{String}, String}, output_path::String, lon_bounds::Tuple{<:Real, <:Real}, lat_bounds::Tuple{<:Real, <:Real})::Nothing
    if typeof(relevant_fields) == String
      dataset = NCDataset(relevant_fields)
    else
      dataset = NCDataset(relevant_fields, aggdim = "")
    end
  
     
    geo_bnds = GeographicBounds(lon_bounds, lat_bounds, dataset)
      
    lon_size = size(geo_bnds.remaining_lon_values, 1)
    lat_size = size(geo_bnds.remaining_lat_values, 1)

    # these variables are used for calculation of pressure levels at each specific lat, lon, time coordinate: p = ap + b * ps
    ap = dataset[:ap][:]
    b = dataset[:b][:]
    time_size = size(dataset[:time], 1)
        

    result_data::Array{Union{Float64, Missing}, 3} = zeros(lon_size, lat_size, time_size)
    println("Time used for calculating the IVT field: ")
    @time begin
      for time in 1:time_size

        hus_data = load_data_in_geo_bounds(dataset, :hus, geo_bnds, :, time)
        ua_data = load_data_in_geo_bounds(dataset, :ua, geo_bnds, :, time)
        va_data = load_data_in_geo_bounds(dataset, :va, geo_bnds, :, time)
        ps_data = load_data_in_geo_bounds(dataset, :ps, geo_bnds, time)

        
        
        Threads.@threads for lat in 1:lat_size
          for lon in 1:lon_size
            
            hus_column::Vector{Union{Float32, Missing}} = hus_data[lon, lat, :]
            ua_column::Vector{Union{Float32, Missing}} = ua_data[lon, lat, :]
            va_column::Vector{Union{Float32, Missing}} = va_data[lon, lat, :]

            ps = ps_data[lon, lat]
            pressure_levels = ap + b * ps
  
            vertical_column_data = VerticalColumnData(hus_column, ua_column, va_column, pressure_levels, ps)
  
            result_data[lon, lat, time] = ivt_of_column(vertical_column_data)
  
          end
        end
      end
      
    end
    
    # copy the data to a new Ordered Dict since a direct pass would only pass the reference to the old atrrib dict -> ERROR!
    new_attrib = OrderedDict()

    for (k, v) in dataset.attrib
      new_attrib[k] = v 
    end
        
      
      # then rewrite some fields
    new_attrib["history"] = "Used an IVT julia script found here: https://github.com/yum-yab/master_thesis/tree/main/preprocessing" * new_attrib["history"]
    new_attrib["title"] = "IVT field of europe and the north-east Atlantic"
  
    NCDataset(output_path, "c", attrib = new_attrib) do ds
  
      defVar(ds, "lon", geo_bnds.remaining_lon_values, ("lon",), attrib = OrderedDict(
        "units" => "degrees east",
        "standard_name" => "longitude"
      ))
      defVar(ds, "lat", geo_bnds.remaining_lat_values, ("lat",), attrib = OrderedDict(
        "units" => "degrees north",
        "standard_name" => "latitude"
      ))
      defVar(ds, "time", dataset[:time][:], ("time",), attrib = OrderedDict(
        "units" => "days since 1850-1-1 00:00:00",
        "standard_name" => "time"
      ))
      defVar(ds, "ivt", result_data, ("lon", "lat", "time"), attrib = OrderedDict(
        "units" => "kg/ms",
        "comments" => "This field is NOT normalized over time or temperature or whatever"
      ))
      
      
    end
  
    close(dataset)
  
    return 
  end
  
  
  
end


module preprocessing

  using NCDatasets
  using Distributed
  using NCDatasets
  using DataStructures
  
  include("IVT.jl")
  using .IVT


  export generate_ivt_field, GeographicBounds, load_data_in_geo_bounds, load_data_single_dataset
  
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
  function load_data_single_dataset{T <: AbstractFloat}(::Type{T}, dataset_path::String, field_id::String, geo_bnds::GeographicBounds, indices...)::Array{T}
  
    lon_normal_range = geo_bnds.lon_bounds[1] < geo_bnds.lon_bounds[2]
    lat_normal_range = geo_bnds.lat_bounds[1] < geo_bnds.lat_bounds[2]
  
    # first load the file information
    dims = Vector{Int}()
  
    NCDataset(dataset_path) do ds
  
      append!(dims, collect(size(ds[field_id])))
      
    end
    
  
    function build_unnormal_ranges(normal_bounds::Tuple{<:Int, <:Int}, unnormal_bounds::Tuple{<:Int, <:Int}, unnormal_index)::Tuple{Vector{UnitRange}, Vector{UnitRange}}
  
      range1 = Vector{UnitRange{Int}}()
      range2 = Vector{UnitRange{Int}}()
  
      other_dims = [typeof(indices[i]) == Colon ? UnitRange(1:dims[2 + i]) : indices[i] for i in eachindex(indices)]
  
      if unnormal_index == 1
        push!(range1, unnormal_bounds[1]:dims[unnormal_index], normal_bounds[1]:normal_bounds[2], other_dims...)
        push!(range2, 1:unnormal_bounds[2], normal_bounds[1]:normal_bounds[2], other_dims...)
      else
        push!(range1, normal_bounds[1]:normal_bounds[2], unnormal_bounds[1]:dims[unnormal_index], other_dims...)
        push!(range2, normal_bounds[1]:normal_bounds[2], 1:unnormal_bounds[2], other_dims...)
      end
  
      return range1, range2
    end
  
    function concat_ranges(range1, range2, dimension)
  
      data1 = zeros(T, length.(range1)...)
      data1 = zeros(T, length.(range2)...)
  
      NCDataset(dataset_path) do 
        NCDatasets.load!(variable(ds, field_id), data1, range1...)
        NCDatasets.load!(variable(ds, field_id), data2, range2...)
      end
      
      return cat(data1, data2; dimension)
    end
  
  
    if lon_normal_range & lat_normal_range
      ranges = [geo_bnds.lon_indices[1]:geo_bnds.lon_indices[2], geo_bnds.lat_indices[1]:geo_bnds.lat_indices[2], indices...]
      data = zeros(T, ranges...)
  
      NCDataset(dataset_path) do ds
        NCDatasets.load!(variable(ds, field_id), data, ranges...)
      end
      return data
  
    elseif !lon_normal_range & lat_normal_range
  
      (range1, range2) = build_unnormal_ranges(geo_bnds.lat_indices, geo_bnds.lon_indices, 1)
      
      return concat_ranges(range1, range2, 1)
    
    elseif lon_normal_range & !lat_normal_range
      
      (range1, range2) = build_unnormal_ranges(geo_bnds.lon_indices, geo_bnds.lat_indices, 2)
      return concat_ranges(range1, range2, 2)
    else
      # last case is both are over 
  
      # first the left side 
  
      (range_left_1, range_left_2) = build_unnormal_ranges(geo_bnds.lon_indices[1]:dims[1], geo_bnds.lat_indices, 2)
      left_side_data = concat_ranges(range_left_1, range_left_2, 2) 
  
      (range_right_1, range_right_2) = build_unnormal_ranges(1:geo_bnds.lon_indices[2], geo_bnds.lat_indices, 2)
      right_side_data = concat_ranges(range_right_1, range_right_2, 2)
  
      return vcat(left_side_data, right_side_data)
    end
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
    println("Time used for loading the data: ") 
    @time begin
      hus_data = load_data_in_geo_bounds(dataset, :hus, geo_bnds, :, :)
      ua_data = load_data_in_geo_bounds(dataset, :ua, geo_bnds, :, :)
      va_data = load_data_in_geo_bounds(dataset, :va, geo_bnds, :, :)
      ps_data = load_data_in_geo_bounds(dataset, :ps, geo_bnds, :)
      
      lon_size = size(hus_data, 1)
      lat_size = size(hus_data, 2)

      # these variables are used for calculation of pressure levels at each specific lat, lon, time coordinate: p = ap + b * ps
      ap = dataset[:ap][:]
      b = dataset[:b][:]
      time_size = size(dataset[:time], 1)
        
    end

    result_data::Array{Union{Float64, Missing}, 3} = zeros(lon_size, lat_size, time_size)
    println("Time used for calculating the IVT field: ")
    @time begin
      Threads.@threads for time in 1:time_size
        for lat in 1:lat_size
          for lon in 1:lon_size
            
            hus_column::Vector{Union{Float32, Missing}} = hus_data[lon, lat, :, time]
            ua_column::Vector{Union{Float32, Missing}} = ua_data[lon, lat, :, time]
            va_column::Vector{Union{Float32, Missing}} = va_data[lon, lat, :, time]

            ps = ps_data[lon, lat, time]
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

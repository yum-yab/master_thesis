module DataLoading
using NCDatasets
using NetCDF

export load_variable_data_in_bounds, GeographicBounds, load_data_in_geo_bounds

struct GeographicBounds
  """This struct sets the geographic boundaries of the datasets. """
  lon_bounds::Tuple{<:Real,<:Real}
  lon_indices::Tuple{<:Int,<:Int}
  remaining_lon_values::Vector{<:Real}
  lat_bounds::Tuple{<:Real,<:Real}
  lat_indices::Tuple{<:Int,<:Int}
  remaining_lat_values::Vector{<:Real}
end

function GeographicBounds(lon_bounds::Tuple{<:Real,<:Real}, lat_bounds::Tuple{<:Real,<:Real}, lon_values::Vector{<:Real}, lat_values::Vector{<:Real})

  function get_overflowing_result(value_bounds, values)::Tuple{Tuple{<:Int,<:Int},Vector{<:Real}}

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

function GeographicBounds(lon_bounds::Tuple{<:Real,<:Real}, lat_bounds::Tuple{<:Real,<:Real}, dataset, lon_id=:lon, lat_id=:lat)
  lon_vals = dataset[lon_id][:]
  lat_vals = dataset[lat_id][:]

  return GeographicBounds(lon_bounds, lat_bounds, lon_vals, lat_vals)
end

function GeographicBounds(lon_bounds::Tuple{<:Real,<:Real}, lat_bounds::Tuple{<:Real,<:Real}, dataset_path::String, lon_id="lon", lat_id="lat")
  lon_vals = ncread(dataset_path, lon_id)
  lat_vals = ncread(dataset_path, lat_id)

  return GeographicBounds(lon_bounds, lat_bounds, lon_vals, lat_vals)
end

"""This function loads data in given geographic bounds. It supports loading values going 'over the end' like the lon rage from 270-40 NOTE: It is expected that the longitude is given in values from 0-360 deg and lat in range from -90:90"""
# function load_variable_data_in_bounds{T,N}(dataset_path::String, field_id::String, geo_bnds::GeographicBounds, indices...)::Array{T,N} where {T<:AbstractFloat,N<:Int}
function load_variable_data_in_bounds(T, dataset_path::String, field_id::String, geo_bnds::GeographicBounds, indices::Union{Integer, UnitRange, StepRange, Colon}...)

  lon_normal_range = geo_bnds.lon_bounds[1] < geo_bnds.lon_bounds[2]
  lat_normal_range = geo_bnds.lat_bounds[1] < geo_bnds.lat_bounds[2]

  # first load the file information
  dims = Vector{Int}()

  NCDataset(dataset_path) do ds

    append!(dims, collect(size(ds[field_id])))

  end


  function build_unnormal_ranges(normal_bounds::Tuple{<:Int,<:Int}, unnormal_bounds::Tuple{<:Int,<:Int}, unnormal_index)::Tuple{Vector{UnitRange},Vector{UnitRange}}

    range1 = Vector{UnitRange{Int}}()
    range2 = Vector{UnitRange{Int}}()

    other_dims = [isa(indices[i], Colon) ? UnitRange(1:dims[2+i]) : indices[i] for i in eachindex(indices)]

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
    data2 = zeros(T, length.(range2)...)

    NCDataset(dataset_path) do ds
      NCDatasets.load!(variable(ds, field_id), data1, range1...)
      NCDatasets.load!(variable(ds, field_id), data2, range2...)
    end

    return cat(data1, data2; dims = dimension)
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
function load_data_in_geo_bounds(dataset, field_id::Union{String,Symbol,Missing}, geo_bnds::GeographicBounds, indices...)::Array

  lon_normal_range = geo_bnds.lon_bounds[1] < geo_bnds.lon_bounds[2]
  lat_normal_range = geo_bnds.lat_bounds[1] < geo_bnds.lat_bounds[2]


  if lon_normal_range & lat_normal_range
    return dataset[field_id][geo_bnds.lon_indices[1]:geo_bnds.lon_indices[2], geo_bnds.lat_indices[1]:geo_bnds.lat_indices[2], indices...]
  elseif !lon_normal_range & lat_normal_range

    lon_first = dataset[field_id][geo_bnds.lon_indices[1]:end, geo_bnds.lat_indices[1]:geo_bnds.lat_indices[2], indices...]
    lon_second = dataset[field_id][1:geo_bnds.lon_indices[2], geo_bnds.lat_indices[1]:geo_bnds.lat_indices[2], indices...]

    return vcat(lon_first, lon_second)
  elseif lon_normal_range & !lat_normal_range

    lat_first = dataset[field_id][geo_bnds.lon_indices[1]:geo_bnds.lon_indices[2], geo_bnds.lat_indices[1]:end, indices...]
    lat_second = dataset[field_id][geo_bnds.lon_indices[1]:geo_bnds.lon_indices[2], 1:geo_bnds.lat_indices[2], indices...]

    return hcat(lat_first, lat_second)
  else
    # last case is both are over 
    lon_f_lat_f = dataset[field_id][geo_bnds.lon_indices[1]:end, geo_bnds.lat_indices[1]:end, indices...]
    lon_f_lat_s = dataset[field_id][geo_bnds.lon_indices[1]:end, 1:geo_bnds.lat_indices[2], indices...]

    lon_s_lat_f = dataset[field_id][1:geo_bnds.lon_indices[2], geo_bnds.lat_indices[1]:end, indices...]
    lon_s_lat_s = dataset[field_id][1:geo_bnds.lon_indices[2], 1:geo_bnds.lat_indices[2], indices...]

    return vcat(hcat(lon_f_lat_f, lon_f_lat_s), hcat(lon_s_lat_f, lon_s_lat_s))
  end
end


end


module Test
using NCDatasets

include("generate_ivt_fields.jl")

using .preprocessing

export load_data_single_dataset

"""This function loads data in given geographic bounds. It supports loading values going 'over the end' like the lon rage from 270-40 NOTE: It is expected that the longitude is given in values from 0-360 deg and lat in range from -90:90""" 
function load_data_single_dataset(dataset_path, field_id::String, geo_bnds, indices...)::Array

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

end


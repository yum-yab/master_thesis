
module preprocessing

  using NCDatasets
  using Distributed
  using NCDatasets
  using DataStructures
  
  include("IVT.jl")
  using .IVT


  export generate_ivt_field, generate_ivt_fields_for_ssp  
  
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
  
    first_lon_id = findfirst(x -> x >= lon_bounds[1], lon_values)
    second_lon_id = findlast(x -> x <= lon_bounds[2], lon_values)

    first_lat_id = findfirst(x -> x >= lat_bounds[1], lat_values)
    second_lat_id = findlast(x -> x <= lat_bounds[2], lat_values)
    
    function get_overflowing_result(bounds, values)
      if bounds[1] < bounds[2]
        remaining_values = values[bounds[1]:bounds[2]]
      else
        remaining_values = vcat(values[bounds[1]:end], values[1:bounds[2]]) 
      end
      return remaining_values
    end
    
    return GeographicBounds(lon_bounds, (first_lon_id, second_lon_id), get_overflowing_result(lon_bounds, lon_values), lat_bounds, (first_lat_id, second_lat_id), get_overflowing_result(lat_bounds, lat_values))
  end

  function GeographicBounds(lon_bounds::Tuple{<:Real, <:Real}, lat_bounds::Tuple{<:Real, <:Real}, dataset, lon_id = :lon, lat_id = :lat)
    lon_vals = dataset[lon_id][:]
    lat_vals = dataset[lat_id][:]

    return GeographicBounds(lon_bounds, lat_bounds, lon_vals, lat_vals)
  end
  

  
  function find_common_versions(member_path::String, field_ids::Vector{String}, time_res_id::String)::Union{String, Nothing}
  
    available_versions = [Set(readdir(joinpath(member_path, time_res_id, fid, "gn"))) for fid in field_ids]
  
    common_versions = Vector(intersect(available_versions...))
    
    
    if isempty(common_versions)
      return nothing 
    else
      # get the latest common version
      return sort(common_versions, rev = true)[1]
    end
  end
  
  function timestamp_from_nc_file(nc_file_name::String)::String
    
    return split(split(nc_file_name, "_")[end], ".")[1]
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
    
    lat_first = dataset[field_id][geo_bnds.lon_bounds[1]:geo_bnds.lon_indices[2], geo_bnds.lat_indices[1]:end,indices...]
    lat_second = dataset[field_id][geo_bnds.lon_bounds[1]:geo_bnds.lon_indices[2], 1:geo_bnds.lat_indices[2],indices...]
    
    return hcat(lat_first, lat_second)
  else
    # last case is both are over 
    lon_f_lat_f = dataset[field_id][geo_bnds.lon_bounds[1]:end, geo_bnds.lat_indices[1]:end,indices...]
    lon_f_lat_s = dataset[field_id][geo_bnds.lon_bounds[1]:end, 1:geo_bnds.lat_indices[2],indices...]
    lon_first = cat(lon_f_lat_f, lon_f_lat_s; dims = 2)

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

    result_data::Array{Union{Float64, Missing}, 3} = zeros(lon_size, lat_size, time_size)
    
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

    new_attrib = OrderedDict()

    for (k, v) in dataset.attrib
      new_attrib[k] = v 
    end
        
      
      # then rewrite some fields
    new_attrib["history"] = "Used an IVT julia script found here: https://github.com/yum-yab/master_thesis/tree/main/preprocessing" * new_attrib["history"]
    new_attrib["title"] = "IVT field of europe and the north-east Atlantic"
  
    NCDataset(output_path, "c", attrib = new_attrib) do ds
  
      # defVar(ds, "lon", geo_bnds.remaining_lon_values, ("lon"), attrib = OrderedDict(
      #   "units" => "degrees east",
      #   "standard_name" => "longitude"
      # ))
      # defVar(ds, "lat", geo_bnds.remaining_lat_values, ("lat"), attrib = OrderedDict(
      #   "units" => "degrees north",
      #   "standard_name" => "latitude"
      # ))
      # defVar(ds, "time", dataset[:time][:], ("time"), attrib = OrderedDict(
      #   "units" => "days since 1850-1-1 00:00:00",
      #   "standard_name" => "time"
      # ))
      defVar(ds, "ivt", result_data, ("lon", "lat", "time"), attrib = OrderedDict(
        "units" => "kg/ms",
        "comments" => "This field is NOT normalized over time or temperature or whatever"
      ))
      
      
    end
  
    close(dataset)
  
    return 
  end
  
  
  
  
  function generate_ivt_fields_for_ssp(base_path::String, ssp_id::String, target_base_path::String, time_res_id::String = "6hrLev", dry_run::Bool = false, cdo_ntreads::Union{Nothing, Int} = nothing)::Nothing
    
    scenario_path = joinpath(base_path, ssp_id)
  
    field_ids = ["hus", "ua", "va"]
  
    member_paths = filter(p -> isdir(p), readdir(scenario_path, join=true))
    
    cdo_arguments::Vector{String} = []
  
    if dry_run
      push!(cdo_arguments, "-T") 
    end
    if !isnothing(cdo_ntreads)
      push!("-P $(cdo_ntreads)")
    end  
  
    for member_path in member_paths
  
      member_id = basename(member_path)
  
      common_version = find_common_versions(member_path, field_ids, time_res_id)
  
      if isnothing(common_version)
        print("No common version found for scenario $ssp_id and member $member_id")
        return nothing
      end
  
      # now here we create the fields 
  
      variable_paths = [joinpath(member_path, time_res_id, variable, "gn", common_version) for variable in field_ids]
  
      all_files_to_merge = collect(zip(nc_files_in_directory.(variable_paths)...))
      
      target_path = joinpath(target_base_path, ssp_id, member_id)
  
      # first generate all merged files in parallel
      pmap(all_files_to_merge) do files_to_merge
  
        timestamp = timestamp_from_nc_file(basename(files_to_merge[1]))
  
        cdo_result_path = joinpath(target_path, "cdo_process_result_$timestamp.nc")
        
        println(prepare_merged_file(files_to_merge, cdo_result_path, cdo_options = cdo_arguments))
        
      end
  
      # now calculate for each timestamped file the ivt field
      
      for merged_file in nc_files_in_directory(target_path)
        
  
  
      end
  
    end
  
  
  
  end
end

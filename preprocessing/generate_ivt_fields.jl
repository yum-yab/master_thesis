
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
    lat_bounds::Tuple{<:Real, <:Real}
  end
  

  function get_range_of_bounds(geo_bnds::GeographicBounds, values::Vector{<:Real}, bound_type::Symbol)
    if bound_type == :lon
      bounds = geo_bnds.lon_bounds
    elseif bound_type == :lat
      bounds = geo_bnds.lat_bounds
    else
      throw(ErrorException("Can't use symbol $bound_type, please use :lon or :lat"))
    end
    
    if bounds[1] < bounds[2]
      first_id = findfirst(x -> x >= bounds[1], values)
      last_id = findlast(x -> x <= bounds[2], values)
      return first_id:last_id
    else
      
      return values[(values .>= bounds[1]) .| (values .<= bounds[2])]  
    end
  end
  

  """
      prepare_merged_file(input_files::Vector{String}, output_file::String)::String
  
  Merges input files, cuts out lon and lat bounds and interpolates the hybrid sigma pressure levels to the configured pressure levels. Result is written to output file.
  """
  function prepare_merged_file(input_files::Vector{String}, output_file::String, pressure_levels::Vector{Int} = round.(Int, LinRange(100000, 0, 47)), cdo_options = [])::String
  
    LON_BOUNDS=(270, 40)
    LAT_BOUNDS=(20, 80)
    
    cdo_command_args = ["cdo", cdo_options..., "-ml2pl,$(join(string.(pressure_levels), ","))", "-sellonlatbox,$(LON_BOUNDS[1]),$(LON_BOUNDS[2]),$(LAT_BOUNDS[1]),$(LAT_BOUNDS[2])", "-selvar,hus,ua,va", "-merge", input_files..., output_file]
    
    cdo_cmd = Cmd(cdo_command_args)
    
    return run(cdo_cmd, wait=true)
  end
  
  function nc_files_in_directory(directory_path::String)::Vector{String}
  
    return filter(path -> endswith(path, ".nc"), readdir(directory_path, join = true))
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
  
  
  
  
  function generate_ivt_field(merged_dataset_path::String, output_path::String)::Nothing
  
    merged_dataset = NCDataset(merged_dataset_path)
  
    lon_size = size(merged_dataset[:lon], 1)
  
    time_size = size(merged_dataset[:time], 1)
    
    lat_size = size(merged_dataset[:lat], 1)
  
    levs = merged_dataset[:lev][:]
    hus_data = merged_dataset[:hus][:, :, :, :]
    ua_data = merged_dataset[:ua][:, :, :, :]
    va_data = merged_dataset[:va][:, :, :, :]
    ps_data = merged_dataset[:ps][:, :, :]
  
    result_data::Array{Union{Float32, Missing}, 3} = zeros(Float32, lon_size, lat_size, time_size)
    
    Threads.@threads for time in 1:time_size
      for lat in 1:lat_size
        for lon in 1:lon_size
          
          hus_column::Vector{Union{Float32, Missing}} = hus_data[lon, lat, :, time]
          ua_column::Vector{Union{Float32, Missing}} = ua_data[lon, lat, :, time]
          va_column::Vector{Union{Float32, Missing}} = va_data[lon, lat, :, time]

          ps::Float32 = ps_data[lon, lat, time]
  
          vertical_column_data = VerticalColumnData(hus_column, ua_column, va_column, ps)
  
          result_data[lon, lat, time] = ivt_of_column(vertical_column_data, levs)
  
        end
      end
    end
    new_attrib = OrderedDict()

    for (k, v) in merged_dataset.attrib
      new_attrib[k] = v 
    end
        
      
      # then rewrite some fields
    new_attrib["history"] = "Used an IVT julia script found here: https://github.com/yum-yab/master_thesis/tree/main/preprocessing" * new_attrib["history"]
    new_attrib["title"] = "IVT field of europe and the north-east Atlantic"
  
    NCDataset(output_path, "c", attrib = new_attrib) do ds
  
      defVar(ds, "ivt", result_data, ("lon", "lat", "time"), attrib = OrderedDict(
        "units" => "kg/ms",
        "comments" => "This field is NOT normalized over time or temperature or whatever"
      ))
      
    end
  
    
  
    close(merged_dataset)
  
    
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

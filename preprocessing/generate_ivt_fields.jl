using NCDatasets
using Distributed
using NCDatasets
using Threads
using DataStructures

include("IVT.jl")
using IVT


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




function generate_ivt_field(merged_fields_path::String, output_path::String)::Nothing

  merged_dataset = NCDataset(merged_fields_path)

  lon_size = size(merged_dataset[:lon], 1)

  time_size = size(merged_dataset[:time], 1)
  
  lat_size = size(merged_dataset[:lat], 1)

  plevs = merged_dataset[:plev][:]
  hus_data = merged_dataset[:hus][:, :, :, :]
  ua_data = merged_dataset[:hus][:, :, :, :]
  va_data = merged_dataset[:hus][:, :, :, :]

  result_data::Array{Union{Float32, Missing}, 3} = zeros(Float32, lon_size, lat_size, time_size)
  
  Threads.@threads for time in 1:time_size
    for lat in 1:lat_size
      for lon in 1:lon_size
        
        hus_column = hus_data[lon, lat, :, time]
        ua_column = ua_data[lon, lat, :, time]
        va_column = va_data[lon, lat, :, time]

        plev_data = map(t -> PressureLevelData(t...), zip(hus_column, ua_column, va_column, plevs))

        result_data[lon, lat, time] = ivt_of_column(plev_data)

      end
    end
  end

  NCDataset(output_path, "c") do ds

    # first copy the attribs as CDO would
    attribs = merged_dataset.attrib 
    
    # then rewrite some fields
    attribs["history"] = "Used an IVT julia script found here: https://github.com/yum-yab/master_thesis/tree/main/preprocessing" * attribs["history"]
    attribs["title"] = "IVT field of europe and the north-east Atlantic"

    ds.attrib = attribs

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

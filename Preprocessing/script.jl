using Distributed
using NCDatasets
addprocs(3, exeflags="--project=$(Base.active_project())")
@everywhere using Preprocessing


function find_common_versions(member_path::String, field_ids::Vector{String}, time_res_id::String)::Union{String, Nothing}

  available_versions = [Set(readdir(joinpath(member_path, time_res_id, fid, "gn"))) for fid in field_ids]

  common_versions = collect(intersect(available_versions...))
  
  
  if isempty(common_versions)
    return nothing 
  else
    # get the latest common version
    return sort(common_versions, rev = true)[1]
  end
end

function nc_files_in_directory(directory_path::String)::Vector{String}

  return filter(path -> endswith(path, ".nc"), readdir(directory_path, join = true, sort = true))
end 

function get_id_to_path_mapping(member_path::String, field_ids::Vector{String}, time_res_id::String)::Union{Vector{Dict{String, String}}, Nothing}
  
  id_to_files = Dict([fid => nc_files_in_directory(readdir(joinpath(member_path, time_res_id, fid, "gn"), join = true)[1]) for fid in field_ids])
  
  # it is assumed no file is missing or this breaks
  all_mappings_for_member = Dict.(zip([[variable => file for file in files] for (variable, files) in id_to_files]...))

  
  return all_mappings_for_member
end

function timestamp_from_nc_file(nc_file_name::String)::String
  
  return split(split(nc_file_name, "_")[end], ".")[1]
end
  



function process_members(process_fun::Function, base_path::String, ssp_id::String, field_ids::Vector{String}; time_res_id::String = "6hrLev", silent::Bool = false)
  scenario_path = joinpath(base_path, ssp_id)

  member_paths = filter(p -> isdir(p), readdir(scenario_path, join=true))
  
  for member_path in member_paths
    member_id = basename(member_path)
    if !silent
      println("Handling member $member_id")
    end
    id_to_path_mappings = get_id_to_path_mapping(member_path, field_ids, time_res_id)

    for id_to_file_mapping in id_to_path_mappings
      process_fun(member_id, id_to_file_mapping)
    end

  end
  
end

function get_id_to_file_mappings(base_path::String, ssp_id::String, field_ids::Vector{String}; time_res_id::String = "6hrLev", silent::Bool = false)::Vector{Dict{String, String}}
  
  result = Vector{Dict{String, String}}()
  process_members(base_path, ssp_id, field_ids; time_res_id = time_res_id, silent = silent) do member_id, id_to_file_mapping
    push!(result, id_to_file_mapping) 
  end
  return result 
end

function generate_ivt_fields_for_ssp(base_path::String, ssp_id::String, target_base_path::String, lon_bnds::Tuple{<:Real, <:Real}, lat_bnds::Tuple{<:Real, <:Real}; time_res_id::String = "6hrLev", overwrite_existing::Bool = true, dry_run::Bool = true)::Nothing
  

  field_ids = ["hus", "ua", "va", "ps"]

  process_members(base_path, ssp_id, field_ids; time_res_id = time_res_id) do member_id, id_to_file_mapping
      timestamp = timestamp_from_nc_file(basename(id_to_file_mapping[field_ids[1]]))
      target_file = joinpath(target_base_path, ssp_id, member_id, "ivt_$(ssp_id)_$(member_id)_$timestamp.nc")

      

      if dry_run
        println("Would have run IVT generation on files: $(id_to_file_mapping) with output being: $target_file")
      else
        if overwrite_existing | !isfile(target_file)
          geo_bounds = DataLoading.GeographicBounds(lon_bnds, lat_bnds, id_to_file_mapping["hus"])
          (data_eastwards, data_northwards, data_norm) = generate_ivt_field(id_to_file_mapping, geo_bounds)
          println("Time it takes saving the data to disk: ")
          @time write_ivt_dataset(id_to_file_mapping[field_ids[1]], geo_bounds, target_file, data_eastwards, data_northwards, data_norm)
        else
          println("Skipped creation of file $target_file: Already existing!") 
        end
      end
      flush(stdout)
    
  end
end

function generate_ivt_fields_xarray(base_path::String, ssp_id::String, target_base_path::String, lon_bnds::Tuple{<:Real, <:Real}, lat_bnds::Tuple{<:Real, <:Real}; time_res_id::String = "6hrLev", overwrite_existing::Bool = true, dry_run::Bool = true)::Nothing
  

  field_ids = ["hus", "ua", "va"]

  process_members(base_path, ssp_id, field_ids; time_res_id = time_res_id) do member_id, id_to_file_mapping
      timestamp = timestamp_from_nc_file(basename(id_to_file_mapping[field_ids[1]]))
      target_file = joinpath(target_base_path, ssp_id, member_id, "ivt_$(ssp_id)_$(member_id)_$timestamp.nc")

      if dry_run
        println("Would have run IVT generation on files: $(id_to_file_mapping) with output being: $target_file")
      else
        if overwrite_existing | !isfile(target_file)
          
          println("Time it took loading the data:")
          
          @time data_dict = XarrayDataLoading.parallel_loading_of_datasets(id_to_file_mapping)
          geo_bounds = DataLoading.GeographicBounds(lon_bnds, lat_bnds, id_to_file_mapping["hus"])
          (data_eastwards, data_northwards, data_norm) = generate_ivt_field(data_dict)
          println("Time it takes saving the data to disk: ")
          @time write_ivt_dataset(id_to_file_mapping[field_ids[1]], geo_bounds, target_file, data_eastwards, data_northwards, data_norm)
        else
          println("Skipped creation of file $target_file: Already existing!") 
        end
      end
      flush(stdout)
    
  end
end


function generate_cdo_preprocessing_commands(base_path::String, ssp_id::String, target_base_path::String, lon_bnds::Tuple{<:Real, <:Real}, lat_bnds::Tuple{<:Real, <:Real}; time_res_id::String = "6hrLev", silent::Bool = false, overwrite_existing::Bool = true, cdo_flags::Vector{String} = Vector{String}())

  field_ids = ["hus", "ua", "va"]
  field_id_selection = join(field_ids, ",")
  
  cdo_commands = Vector{Cmd}()

  process_members(base_path, ssp_id, field_ids; time_res_id = time_res_id, silent = silent) do member_id, id_to_file_mapping 
    
    timestamp = timestamp_from_nc_file(basename(id_to_file_mapping[field_ids[1]]))

    target_file = joinpath(target_base_path, ssp_id, member_id, "merged-preprocessed-data_$(ssp_id)_$(member_id)_$timestamp.nc")
    mkpath(dirname(target_file))

    source_files = [filepath for (_, filepath) in id_to_file_mapping]
    
    cdo_cmd = `cdo $cdo_flags -selvar,$field_id_selection -sellonlatbox,$(lon_bnds[1]),$(lon_bnds[2]),$(lat_bnds[1]),$(lat_bnds[2]) -merge $source_files $target_file`

    if !overwrite_existing || !isfile(target_file)
      push!(cdo_commands, cdo_cmd)
    end
  
  end

  return cdo_commands
end
function generate_cdo_preprocessing_commands(base_path::String, ssp_id::String, target_base_path::String, lon_bnds::Tuple{<:Real, <:Real}, lat_bnds::Tuple{<:Real, <:Real}; time_res_id::String = "6hrLev", silent::Bool = false, overwrite_existing::Bool = true, cdo_flags::Vector{String} = Vector{String}())

  field_ids = ["hus", "ua", "va"]
  field_id_selection = join(field_ids, ",")
  
  cdo_commands = Vector{Cmd}()

  process_members(base_path, ssp_id, field_ids; time_res_id = time_res_id, silent = silent) do member_id, id_to_file_mapping 
    
    timestamp = timestamp_from_nc_file(basename(id_to_file_mapping[field_ids[1]]))

    target_file = joinpath(target_base_path, ssp_id, member_id, "merged-preprocessed-data_$(ssp_id)_$(member_id)_$timestamp.nc")
    mkpath(dirname(target_file))

    source_files = [filepath for (_, filepath) in id_to_file_mapping]
    
    cdo_cmd = `cdo $cdo_flags -selvar,$field_id_selection -sellonlatbox,$(lon_bnds[1]),$(lon_bnds[2]),$(lat_bnds[1]),$(lat_bnds[2]) -merge $source_files $target_file`

    if !overwrite_existing || !isfile(target_file)
      push!(cdo_commands, cdo_cmd)
    end
  
  end

  return cdo_commands
end


function generate_ivt_from_preprocessed_files(base_path::String, ssp_id::String, target_base_path::String; overwrite_existing::Bool = true)

  member_paths = filter(p -> isdir(p), readdir(joinpath(base_path, ssp_id), join=true))
  field_ids = ["hus", "ua", "va", "ps"]
  for member_path in member_paths

    member_id = basename(member_path)

    merge_member_files = filter(p -> isfile(p), readdir(member_path, join=true))
    
    for merged_file in merge_member_files

      id_to_path_dict = Dict([fid => merged_file for fid in field_ids])
      timestamp = timestamp_from_nc_file(basename(merged_file))
      target_file = joinpath(target_base_path, ssp_id, member_id, "ivt_$(ssp_id)_$(member_id)_$timestamp.nc")

      

      (data_east, data_north, data_norm) = generate_ivt_field(id_to_path_dict)

      println("Time it takes saving the data to disk: ")
      @time write_ivt_dataset(merged_file, target_file, data_east, data_north, data_norm)

      
    end      
  end
  
end


function main(cfg::Dict{String, Any})
  println("Number of threads available: $(Threads.nthreads())")  
  scenario_ssps = cfg["source_selection"]["ssps"]
  time_res_id = cfg["source_selection"]["time_resolution_id"]
  lon_bounds = Tuple(cfg["dataset"]["lon_bounds"])
  lat_bounds = Tuple(cfg["dataset"]["lat_bounds"])
  scenrio_base_path = cfg["process"]["scenario_path"]
  target_base_path = cfg["process"]["output_path"]
  overwrite_existing = cfg["process"]["overwrite_existing"]

  dry_run = cfg["process"]["dry_run"]

  for ssp in scenario_ssps
    id_to_file_mappings = get_id_to_file_mappings(scenrio_base_path, ssp, ["hus", "ua", "va"]; silent = true)
    target_file = joinpath(target_base_path, "test_xarray_full_script.nc") 
    id_to_file_mapping = id_to_file_mappings[1]
    full_mapping_dict = merge(id_to_file_mapping, Dict("ap" => id_to_file_mapping["hus"], "b" => id_to_file_mapping["hus"]))
    
    println("Time for loading data:")
    @time data_dict = XarrayDataLoading.parallel_loading_of_datasets(full_mapping_dict)
    rmprocs(workers())
    NCDataset(id_to_file_mapping["hus"]) do ds

      data_dict["ps"] = ds["ps"][:, :, :]
      data_dict["ap"] = ds["ap"][:]
      data_dict["b"] = ds["b"][:]
      
      
    end
    geo_bounds = DataLoading.GeographicBounds(lon_bounds, lat_bounds, id_to_file_mapping["hus"])
    (data_eastwards, data_northwards, data_norm) = generate_ivt_field(data_dict)
    println("Time it takes saving the data to disk: ")
    @time write_ivt_dataset(id_to_file_mapping["hus"], geo_bounds, target_file, data_eastwards, data_northwards, data_norm)
 
    # commands = generate_cdo_preprocessing_commands(scenrio_base_path, ssp, target_base_path, lon_bounds, lat_bounds; time_res_id = time_res_id, overwrite_existing = overwrite_existing, silent = true)

    # for cmd in commands
    #     println(replace(string(cmd), "`" => ""))
    # end

  end
  
end

using TOML 

main(TOML.parsefile("./cfg.toml"))

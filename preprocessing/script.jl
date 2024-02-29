include("generate_ivt_fields.jl")

using .preprocessing



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

function timestamp_from_nc_file(nc_file_name::String)::String
  
  return split(split(nc_file_name, "_")[end], ".")[1]
end
  

function nc_files_in_directory(directory_path::String)::Vector{String}

  return filter(path -> endswith(path, ".nc"), readdir(directory_path, join = true))
end 


function generate_ivt_fields_for_ssp(base_path::String, ssp_id::String, target_base_path::String, lon_bnds::Tuple{<:Real, <:Real}, lat_bnds::Tuple{<:Real, <:Real}; time_res_id::String = "6hrLev", overwrite_existing::Bool = true, dry_run::Bool = true)::Nothing
  
  scenario_path = joinpath(base_path, ssp_id)

  field_ids = ["hus", "ua", "va", "ps"]

  member_paths = filter(p -> isdir(p), readdir(scenario_path, join=true))

  for member_path in member_paths
  
    member_id = basename(member_path)
    println("Running IVT field generatioIt is possible to redirect job output to somewhere other than the default location with the --error and --output directives:n for member $member_id")

    common_version = find_common_versions(member_path, field_ids, time_res_id)

    if isnothing(common_version)
      print("No common version found for scenario $ssp_id and member $member_id")
      return nothing
    end

    # now here we create the fields 

    variable_to_member_files = Dict(variable => nc_files_in_directory.(joinpath(member_path, time_res_id, variable, "gn", common_version)) for variable in field_ids)

    # it is assumed no file is missing or this breaks
    all_mappings_for_member = Dict.(zip([[variable => file for file in files] for (variable, files) in variable_to_member_files]...))
    
    target_path = joinpath(target_base_path, ssp_id, member_id)
    mkpath(target_path)
    
    for id_to_file_mapping in all_mappings_for_member
      timestamp = timestamp_from_nc_file(basename(id_to_file_mapping[field_ids[1]]))
      target_file = joinpath(target_path, "ivt_$(ssp_id)_$(member_id)_$timestamp.nc")

      if dry_run
        println("Would have run IVT generation on files: $(id_to_file_mapping) with output being: $target_file")
      else
        if overwrite_existing | !isfile(target_file)
          geo_bounds = preprocessing.DataLoading.GeographicBounds(lon_bnds, lat_bnds, id_to_file_mapping["hus"])
          println("Time it took for the whole ivt field generation of memeber $member_id and timeslice $timestamp :")
          data = generate_ivt_field(id_to_file_mapping, geo_bounds)
          println("Time it takes saving the data to disk: ")
          @time write_ivt_dataset(id_to_file_mapping[field_ids[1]], geo_bounds, target_file, data)
        else
          println("Skipped creation of file $target_file: Already existing!") 
        end
      end
      flush(stdout)
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
    generate_ivt_fields_for_ssp(scenrio_base_path, ssp, target_base_path, lon_bounds, lat_bounds; time_res_id = time_res_id, overwrite_existing = overwrite_existing, dry_run = dry_run)
  end
  
end

using TOML 

main(TOML.parsefile("./cfg.toml"))

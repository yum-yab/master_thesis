


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


function register_sbatch_commands(scenario_base_path, target_base_path, logging_base_path, scenarios...; really_register = true)
  if isempty(PROGRAM_FILE)
    working_dir="."
  else
    working_dir = dirname(abspath(PROGRAM_FILE))
  end
  for scenario in scenarios

    all_members = filter(path -> isdir(path), readdir(joinpath(scenario_base_path, scenario), join = true, sort = true))

    for member_path in all_members
      member_id = basename(member_path)
      jobname = "ivt_gen_$(scenario)_$(member_id)"
      output_log = "$logging_base_path/output_logs/output_$(scenario)_$(member_id).log"
      error_log = "$logging_base_path/error_logs/error_$(scenario)_$(member_id).log"
      sbatch_command=`sbatch --job-name=$jobname --chdir=$working_dir --output=$output_log --error=$error_log run_pythonscript.sh $member_path $target_base_path`
      println(sbatch_command)
      if really_register
      	run(sbatch_command)
      end
    end
  end
  
end




function main()
  
  target_base = "/scratch/b/b382641/ivt_fields_v1"
  scenario_base = "/pool/data/CMIP6/data/ScenarioMIP/MPI-M/MPI-ESM1-2-LR"
  logging_base_path = "/home/b/b382641/workspace/master_thesis/Preprocessing"

  scenarios = ["ssp585"]

  register_sbatch_commands(scenario_base, target_base, logging_base_path, scenarios...; really_register = false)
end

main()

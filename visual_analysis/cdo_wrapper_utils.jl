using Distributed


function run_cdo_command_on_all(base_path, result_basepath, command; run_command = false, overwrite = true)

    for scenario_path in readdir(base_path, join = true, sort = true)
        
        scenario_id = basename(scenario_path)
        
        for member_path in readdir(scenario_path, join = true, sort = true)

            member_id = basename(member_path)
            
            for dataset_path in readdir(member_path, join = true, sort = true)
                
                dataset_name = basename(dataset_path)

                println("Handling member $(member_id) from scenario $(scenario_id)...")

                output_path = joinpath(result_basepath, scenario_id, member_id)

                mkpath(output_path)

                output_file = joinpath(output_path, dataset_name)

                if isfile(output_file) && !overwrite
                    continue
                end


                cdo_cmd = `cdo $command $(dataset_path) $(output_file)`
                
                println(cdo_cmd)
                if run_command
                    run(cdo_cmd)
                end
            end
        end
    end
    
end

function run_cdo_command_on_all_multiprocessing(base_path, result_basepath, command; run_command = false, overwrite = false)

    for scenario_path in readdir(base_path, join = true, sort = true)
        
        scenario_id = basename(scenario_path)

        source_paths = String[]

        output_paths = String[]

        
        for member_path in readdir(scenario_path, join = true, sort = true)

            member_id = basename(member_path)
            
            for dataset_path in readdir(member_path, join = true, sort = true)
                
                dataset_name = basename(dataset_path)

                println("Handling member $(member_id) from scenario $(scenario_id)...")

                output_path = joinpath(result_basepath, scenario_id, member_id)

                mkpath(output_path)

                output_file = joinpath(output_path, dataset_name)
                if isfile(output_file) && !overwrite
                    continue
                end



                push!(source_paths, dataset_path)
                push!(output_paths, output_file)
                
                
                
            end
        end

        pmap([command for _ in eachindex(source_paths)], source_paths, output_paths, [run_command for _ in eachindex(source_paths)]) do command, dataset_path, output_file, run_command
            cdo_cmd = `cdo $command $(dataset_path) $(output_file)`

                println(cdo_cmd, " handled in thread $(myid())")
                if run_command
                    run(cdo_cmd)
                end
        end


    end
    
end

if length(ARGS) == 3
    command = ARGS[1]
    source_path = ARGS[2]
    target_path = ARGS[3]
    run_cdo_command_on_all(source_path, target_path, [command]; run_command = true, overwrite = false)
elseif length(ARGS) == 4
    command = ARGS[1]
    source_path = ARGS[2]
    target_path = ARGS[3]
    num_workers = parse(Int, ARGS[4])

    addprocs(num_workers)

    run_cdo_command_on_all_multiprocessing(source_path, target_path, [command]; run_command = true, overwrite = false)
end
using Distributed


function run_cdo_command_on_all(base_path, result_basepath, command; run_command = false)

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


                cdo_cmd = `cdo $command $(dataset_path) $(output_file)`
                
                println(cdo_cmd)
                if run_command
                    run(cdo_cmd)
                end
            end
        end
    end
    
end

function run_cdo_command_on_all_multiprocessing(base_path, result_basepath, command; run_command = false)

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

if length(ARGS) == 0
    run_cdo_command_on_all("/home/denis/workspace/data/ps_data", "/home/denis/workspace/data/ps_data_daily", ["-daymean"]; run_command = true)
else 
    num_workers = parse(Int, ARGS[1])

    addprocs(num_workers)

    run_cdo_command_on_all_multiprocessing("/home/denis/workspace/data/ps_data", "/home/denis/workspace/data/ps_data_daily_multi", ["-daymean"]; run_command = true)
end
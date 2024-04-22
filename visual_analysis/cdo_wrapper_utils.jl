
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


run_cdo_command_on_all("/home/denis/workspace/data/ivt_fields_v1", "/home/denis/workspace/data/ivt_monthly_mean", ["-monmean"]; run_command = true)
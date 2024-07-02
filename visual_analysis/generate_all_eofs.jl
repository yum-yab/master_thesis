include("utils.jl")

using Dates
using JLD2


function get_sliding_window(time_data, n)

    time_length = length(time_data)

    last_starting_point = time_length - n + 1

    return [i:(i+n-1) for i in 1:last_starting_point]

end

function filter_winter_season(time_element)
    winter_months = [12, 1, 2, 3]
    for wm in winter_months
        if month(time_element) == wm
            return true
        end
    end
    return false
end

function winter_timelimit(timeelement)
    return filter_winter_season(timeelement) && timeelement < DateTime(2100, 12, 31)
end



function generate_access_dict(data_base_path, data_field_id)



    (historical, ssp126, ssp585) = build_ensemble_data(data_base_path, "historical", "ssp126", "ssp585"; file_range_selection=:, data_field_id=data_field_id, member_range=1:50, filterfun=winter_timelimit)


    (piControl,) = build_ensemble_data(data_base_path, "piControl"; file_range_selection=1:13, data_field_id=data_field_id, member_range=1:1, filterfun=winter_timelimit)

    ssp126_full_timeline = concat_ensemble_data(historical, ssp126; id="$(uppercase(data_field_id)) SSP126 full timeline")
    ssp585_full_timeline = concat_ensemble_data(historical, ssp585; id="$(uppercase(data_field_id)) SSP585 full timeline")



    access_dict = Dict(
        "ssp126" => ssp126_full_timeline,
        "ssp585" => ssp585_full_timeline,
        "piControl" => piControl,
        "time" => ssp585_full_timeline.time
    )

    return access_dict
end


function generate_eof_data(
    access_dict,
    variable_label,
    target_base_path;
    nmodes=5,
    season_timescopes=Dict(i => "$(i)seasons" for i in [30, 50]),
    scale_by_sqrt=Dict(false => "nosqrtscale", true => "sqrtscale"),
    scenario_ids=["ssp585", "ssp126"],
    multithreading=false
)

    for (sqrtscale, sqrtlabel) in scale_by_sqrt
        for (scope_size, scope_label) in season_timescopes

            scopes = get_sliding_time_scopes_by_threshold(access_dict["time"], scope_size)

            piControl_eof = calculate_eofs_of_ensemble_fast(
                access_dict["piControl"],
                scopes,
                nmodes;
                center=true,
                align_eofs_with_mean=true,
                norm_withsqrt_timedim=sqrtscale,
                geoweights=true,
                scale_mode=noscaling,
                multithreading=multithreading
            )
            for scenario in scenario_ids

                data_eof = calculate_eofs_of_ensemble_fast(
                    access_dict[scenario],
                    scopes,
                    nmodes;
                    center=true,
                    align_eofs_with_mean=true,
                    norm_withsqrt_timedim=sqrtscale,
                    geoweights=true,
                    scale_mode=noscaling,
                    multithreading=multithreading
                )

                file_name = "$(variable_label)_eofs_$(nmodes)modes_$(scenario)_$(scope_label)_$(sqrtlabel).jld2"

                persisting_dir = joinpath(target_base_path, variable_label, scenario, scope_label)

                mkpath(persisting_dir)


                jldsave(joinpath(persisting_dir, file_name), variable_id=variable_label, scopes=scopes, time=access_dict["time"], piControl=piControl_eof, scenario_data=data_eof)

            end
        end
    end

end

function main(data_base_path, target_base_path)

    # for (variable_id, dirname) in [("pr", "preciptation_data_monthly"), ("ivt", "ivt_fields_monthly"), ("ps", "ps_data_monthly")]
    for (variable_id, dirname) in [("psl", "psl_data_monthly")]

        access_dict = generate_access_dict(joinpath(data_base_path, dirname), variable_id)
        generate_eof_data(
            access_dict,
            variable_id,
            target_base_path;
            nmodes=5,
            season_timescopes=Dict(i => "$(i)seasons" for i in [30, 50]),
            scale_by_sqrt=Dict(true => "sqrtscale"),
            scenario_ids=["ssp585", "ssp126"],
            multithreading=false
        )
    end

end

main(ARGS[1], ARGS[2])

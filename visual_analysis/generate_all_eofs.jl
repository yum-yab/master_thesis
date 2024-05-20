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

function piControl_winter_timelimit(timeelement)
    return filter_winter_season(timeelement) && timeelement < DateTime(2101, 2, 1)
end



function generate_access_dict()
    ps_data_monthly_path = "/mnt/bigdrive/Datasets/master_thesis_data/ps_data_monthly"
    ivt_data_monthly_path = "/mnt/bigdrive/Datasets/master_thesis_data/ivt_fields_monthly"

    (ivt_historical_monthly, ivt_ssp126_monthly, ivt_ssp585_monthly) = build_ensemble_data(ivt_data_monthly_path, "historical", "ssp126", "ssp585"; file_range_selection=:, data_field_id="ivt", member_range=1:50, filterfun=filter_winter_season)

    (ps_historical_monthly, ps_ssp126_monthly, ps_ssp585_monthly) = build_ensemble_data(ps_data_monthly_path, "historical", "ssp126", "ssp585"; file_range_selection=:, data_field_id="ps", member_range=1:50, filterfun=filter_winter_season)

    (ps_piControl,) = build_ensemble_data(ps_data_monthly_path, "piControl"; file_range_selection=1:13, data_field_id="ps", member_range=1:1, filterfun=piControl_winter_timelimit)
    (ivt_piControl,) = build_ensemble_data(ivt_data_monthly_path, "piControl"; file_range_selection=1:13, data_field_id="ivt", member_range=1:1, filterfun=piControl_winter_timelimit)

    ivt_ssp126_full_timeline = concat_ensemble_data(ivt_historical_monthly, ivt_ssp126_monthly; id="IVT SSP126 full timeline")
    ivt_ssp585_full_timeline = concat_ensemble_data(ivt_historical_monthly, ivt_ssp585_monthly; id="IVT SSP585 full timeline")

    ps_ssp126_full_timeline = concat_ensemble_data(ps_historical_monthly, ps_ssp126_monthly; id="PS SSP126 full timeline")
    ps_ssp585_full_timeline = concat_ensemble_data(ps_historical_monthly, ps_ssp585_monthly; id="PS SSP585 full timeline")


    access_dict = Dict(
        "ssp126" => (ivt_ssp126_full_timeline, ps_ssp126_full_timeline),
        "ssp585" => (ivt_ssp585_full_timeline, ps_ssp585_full_timeline),
        "ps_piControl" => ps_piControl,
        "ivt_piControl" => ivt_piControl
    )

    return access_dict, ivt_ssp585_full_timeline.time
end




function main()

    access_dict, time_data = generate_access_dict()
    scale_by_sqrt = Dict(
      # false => "nosqrtscale", 
      true => "sqrtscale"
    )

    season_timescopes = Dict(i => "$(i)seasons" for i in [4, 7, 12, 20, 30, 50, 70])

    # month_timescopes = Dict(i => "$(i)months" for i in [8, 15, 30, 50])


    target_dir = "/mnt/bigdrive/Datasets/master_thesis_data/eof_data"

    for (sqrtscale, sqrtlabel) in scale_by_sqrt

        all_scopes = [(season_label, get_sliding_time_scopes_by_threshold(time_data, season_scope)) for (season_scope, season_label) in season_timescopes]

        # append!(all_scopes, [(label, get_sliding_window(time_data, scope)) for (scope, label) in month_timescopes])


        for (scope_label, scopes) in all_scopes

            println("Started eof generation of $sqrtlabel $scope_label")

            ivt_piControl_eof = calculate_eofs_of_ensemble_fast(
                access_dict["ivt_piControl"],
                scopes,
                2;
                center=true,
                align_eofs_with_mean=true,
                norm_withsqrt_timedim=false,
                geoweights=true,
                scale_mode=:singularvals,
            )

            ps_piControl_eof = calculate_eofs_of_ensemble_fast(
                access_dict["ps_piControl"],
                scopes,
                2;
                center=true,
                align_eofs_with_mean=true,
                norm_withsqrt_timedim=false,
                geoweights=true,
                scale_mode=:singularvals,
            )


            for scenario in ["ssp126", "ssp585"]

                ivt_data, ps_data = access_dict[scenario]

                ivt_eof = calculate_eofs_of_ensemble_fast(
                    ivt_data,
                    scopes,
                    2;
                    center=true,
                    align_eofs_with_mean=true,
                    norm_withsqrt_timedim=sqrtscale,
                    geoweights=true,
                    scale_mode=:singularvals,
                )

                ps_eof = calculate_eofs_of_ensemble_fast(
                    ps_data,
                    scopes,
                    2;
                    center=true,
                    align_eofs_with_mean=true,
                    norm_withsqrt_timedim=sqrtscale,
                    geoweights=true,
                    scale_mode=:singularvals,
                )

                file_name = "eofs_$(scenario)_$(scope_label)_$(sqrtlabel).jld2"

                persisting_dir = joinpath(target_dir, scope_label)

                mkpath(persisting_dir)

                jldsave(joinpath(persisting_dir, file_name), scopes=scopes, ivt_piControl=ivt_piControl_eof, ps_piControl=ps_piControl_eof, ivt_eof=ivt_eof, ps_eof=ps_eof)

            end

        end

    end


end

main()

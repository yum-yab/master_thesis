using EmpiricalOrthogonalFunctions
using NCDatasets
using NetCDF
using Dates
using BenchmarkTools
using Statistics
using JLD2
using ProgressBars
using Contour

include("eof.jl")


struct ScenarioData
    name::String
    data::AbstractArray{Union{Missing,<:AbstractFloat},3}
end

struct TimelineData
    lons::Vector{Union{Missing,<:AbstractFloat}}
    lats::Vector{Union{Missing,<:AbstractFloat}}
    time::Vector{Union{Missing,Dates.DateTime}}
    datasets::Vector{ScenarioData}
end


struct EnsembleMember
    id::String
    data::AbstractArray{Union{Missing,<:AbstractFloat},3}
end


struct EnsembleSimulation
    id::String
    lons::Vector{Union{Missing,<:AbstractFloat}}
    lats::Vector{Union{Missing,<:AbstractFloat}}
    time::Vector{Union{Missing,Dates.DateTime}}
    members::Vector{EnsembleMember}
end


struct EOFEnsembleResult
    id::String
    scopes::Vector{UnitRange{Int}}
    ivt_piControl::Vector{EOFResult}
    ps_piControl::Vector{EOFResult}
    ivt_ensemble_eofs::Dict{String,Vector{EOFResult}}
    ps_ensemble_eofs::Dict{String,Vector{EOFResult}}
end

struct EOFEnsemble
    variable_id::String
    time::Vector{DateTime}
    lon::Vector{Union{Missing,<:AbstractFloat}}
    lat::Vector{Union{Missing,<:AbstractFloat}}
    scopes::Vector{UnitRange{Int}}
    piControl::Vector{EOFResult}
    ensemble::Dict{String,Vector{EOFResult}}
end

function get_member_id_string(member_nr::Int)::String
    return "r$(member_nr)i1p1f1"
end

function get_files_of_member(data_path, scenario_id, member_nr)
    return readdir(joinpath(data_path, scenario_id, get_member_id_string(member_nr)), join=true)
end



function quick_unit_lookup(field_id)

    unit_dict = Dict(
        "psl" => "hPa",
        "pr" => "mm/month",
        "ivt" => "kg s-1 m-1"
    )

    return unit_dict[field_id]

end

function get_data(data_path, scenario_id, member_nr; file_range_selection=:, field_id="ivt", unit_scale_factor=1)

    file_paths = get_files_of_member(data_path, scenario_id, member_nr)

    ivt_data = Array{Float64,3}[]



    for file_path in file_paths[file_range_selection]
        ivt_chunk = ncread(file_path, field_id)
        push!(ivt_data, ivt_chunk * unit_scale_factor)
    end

    return cat(ivt_data..., dims=3)
end

function get_time_data(data_path, scenario_id, member_nr; file_range_selection=:)

    file_paths = get_files_of_member(data_path, scenario_id, member_nr)

    time_data = DateTime[]

    for file_path in file_paths[file_range_selection]
        time_chunk = NCDataset(file_path) do ds
            data = ds[:time][:]
            return convert(Vector{DateTime}, data)
        end
        append!(time_data, time_chunk)
    end

    return time_data
end

function get_field(path, field_id)

    data = ncread(path, field_id)

    return data
end

function build_timeline_data(base_path, member, scenarios...; file_range_selection=:, data_field_id="ivt")

    lons = get_field(get_files_of_member(base_path, scenarios[1], member)[1], "lon")
    lats = get_field(get_files_of_member(base_path, scenarios[1], member)[1], "lat")

    time = get_time_data(base_path, scenarios[1], member; file_range_selection=file_range_selection)

    scenarios = map(scenarios) do scenario_id

        data = get_data(base_path, scenario_id, member; file_range_selection=file_range_selection, field_id=data_field_id)

        return ScenarioData(scenario_id, data)
    end

    return TimelineData(lons, lats, time, collect(scenarios))
end

function is_different_month(date1, date2)

    return !(year(date1) == year(date2) && month(date1) == month(date2))
end


# function check_time_alignment(time_axis...)

#     size_set = Set([length(a) for a in time_axis])

#     if length(size_set) != 1
#         return false
#     end

#     for i in range(1, collect(size_set)[1])

function build_ensemble_data_common_scaling(base_path, scenarios...; file_range_selection=:, data_field_id="ivt", member_range=1:50, silent=false, filterfun=nothing)

    (_, scaling, _, _) = get_correct_var_display(data_field_id)
    return build_ensemble_data(base_path, scenarios...; file_range_selection=file_range_selection, data_field_id=data_field_id, member_range=member_range, silent=silent, filterfun=filterfun, unit_scale_factor=scaling)
end


function build_ensemble_data(base_path, scenarios...; file_range_selection=:, data_field_id="ivt", member_range=1:50, silent=false, filterfun=nothing, unit_scale_factor=1)

    lons = ncread(get_files_of_member(base_path, scenarios[1], 1)[1], "lon")
    lats = get_field(get_files_of_member(base_path, scenarios[1], 1)[1], "lat")

    result = EnsembleSimulation[]




    for scenario in scenarios
        if !silent
            println("Handling scenario $scenario ...")
        end

        time = get_time_data(base_path, scenario, 1; file_range_selection=file_range_selection)

        last_value_vec = isnothing(filterfun) || filterfun(time[end]) ? [eachindex(time)[end]] : []

        time_selector = isnothing(filterfun) ? vcat([i for i in eachindex(time)[1:end-1] if is_different_month(time[i], time[i+1])], last_value_vec) : vcat([i for i in eachindex(time)[1:end-1] if filterfun(time[i]) && is_different_month(time[i], time[i+1])], last_value_vec)
        time = time[time_selector]


        ensemble_members = map(member_range) do member_nr

            member_id = get_member_id_string(member_nr)

            data = get_data(base_path, scenario, member_nr; file_range_selection=file_range_selection, field_id=data_field_id, unit_scale_factor=unit_scale_factor)
            return EnsembleMember(member_id, data[:, :, time_selector])
        end

        push!(result, EnsembleSimulation(scenario, lons, lats, time, collect(ensemble_members)))
        flush(stdout)
    end

    return result
end

function concat_ensemble_data(ensemble_simulations::EnsembleSimulation...; id::String=nothing)::EnsembleSimulation

    member_length = Set([length(es.members) for es in ensemble_simulations])

    if isnothing(id)
        id = ensemble_simulations[1].id
    end

    if length(member_length) != 1
        ArgumentError("Simulations have not aligning members: Sizes: $(member_length)")
    end

    first_timestamps = [es.time[1] for es in ensemble_simulations]

    if !issorted(first_timestamps)
        ArgumentError("Simulations have not aligning timescopes: Sizes: $(first_timestamps)")
    end


    new_lats = ensemble_simulations[1].lats
    new_lons = ensemble_simulations[1].lons

    # check overlapping monthly data 

    function get_last_indices_array(ensemble_simulations)

        result = []

        for i in eachindex(ensemble_simulations)

            if i == 1
                continue
            end

            last_time_array = ensemble_simulations[i-1].time
            first_new_date = ensemble_simulations[i].time[1]

            last_index = findlast(last_time_array) do dtime

                return month(dtime) != month(first_new_date)
            end
            push!(result, last_index)
        end

        push!(result, lastindex(ensemble_simulations[end].time))

        return result
    end


    last_indeces_ensemble = get_last_indices_array(ensemble_simulations)

    new_time = cat([es.time[1:last_indeces_ensemble[i]] for (i, es) in enumerate(ensemble_simulations)]..., dims=1)

    members = Vector{EnsembleMember}(undef, length(ensemble_simulations[1].members))

    for index in eachindex(ensemble_simulations[1].members)

        member_id_set = Set([es.members[index].id for es in ensemble_simulations])

        if length(member_id_set) != 1
            ArgumentError("Not aligning members: $member_id_set at index $index")
        end

        members[index] = EnsembleMember(pop!(member_id_set), cat([es.members[index].data[:, :, 1:last_indeces_ensemble[i]] for (i, es) in enumerate(ensemble_simulations)]..., dims=3))

    end


    return EnsembleSimulation(id, new_lons, new_lats, new_time, members)
end

function filter_by_date(fun, timeline_data::TimelineData)::TimelineData

    time_indices = [i for i in eachindex(timeline_data.time) if fun(timeline_data.time[i])]

    transformed_scenarios = [ScenarioData(ds.name, ds.data[:, :, time_indices]) for ds in timeline_data.datasets]

    return TimelineData(timeline_data.lons, timeline_data.lats, timeline_data.time[time_indices], transformed_scenarios)
end

function filter_by_date(fun, ensemble_simulation::EnsembleSimulation)::EnsembleSimulation

    time_indices = [i for i in eachindex(ensemble_simulation.time) if fun(ensemble_simulation.time[i])]

    return EnsembleSimulation(ensemble_simulation.id, ensemble_simulation.lons, ensemble_simulation.lats, ensemble_simulation.time[time_indices], collect(map(ensemble_simulation.members, member -> EnsembleMember(member.id, member.data[:, :, time_indices]))))
end


function align_with_field!(field, alignment_field; dims=1)

    for slice in eachslice(field, dims=dims)
        scalar_product = sum(slice .* alignment_field)
        if scalar_product < 0
            slice *= -1
        end
    end
end



function get_sliding_time_scopes_by_threshold(time_data, n, threshold=Day(100))

    all_starting_points = [1]

    for i in 1:length(time_data)-1
        if time_data[i+1] - time_data[i] > threshold
            push!(all_starting_points, i + 1)
        end
    end

    function get_next_stop(start)

        counter = 1

        for i in start:length(time_data)-1
            if time_data[i+1] - time_data[i] > threshold
                if counter == n
                    return i
                else
                    counter += 1
                end
            end
        end

        return
    end



    return [sp:get_next_stop(sp) for sp in all_starting_points if !isnothing(get_next_stop(sp))]
end

function centralize_tl_data(tl_data; dims=3)

    normalized_ds = ScenarioData[]
    for ds in tl_data.datasets
        res = ds.data .- mean(ds.data, dims=dims)
        push!(normalized_ds, ScenarioData(ds.name, res))
    end

    return TimelineData(tl_data.lons, tl_data.lats, tl_data.time, normalized_ds)
end


function normalize_tl_data(tl_data)

    normalized_ds = ScenarioData[]
    for ds in tl_data.datasets
        min, max = extrema(ds.data)
        res = (ds.data .- min) ./ (max .- min)
        push!(normalized_ds, ScenarioData(ds.name, res))
    end

    return TimelineData(tl_data.lons, tl_data.lats, tl_data.time, normalized_ds)
end



function calculate_eofs_of_ensemble_fast(
    ensemble::EnsembleSimulation,
    chunking,
    nmodes;
    center=true,
    align_eofs_with_mean=true,
    norm_withsqrt_timedim=false,
    geoweights=true,
    scale_mode=noscaling,
    saving_filepath=nothing,
    multithreading=false
)::Dict{String,Vector{EOFResult}}

    result = Dict{String,Vector{EOFResult}}()

    if geoweights
        weights = sqrt.(cos.(deg2rad.(ensemble.lats)))
    else
        weights = nothing
    end

    iter = ProgressBar(ensemble.members)

    for member in iter
        # update description of progress bar
        set_description(iter, "$(ensemble.id) $(member.id)")


        eofs = Vector{EOFResult}(undef, length(chunking))

        if multithreading
            Threads.@threads for idx in eachindex(chunking)
                if scale_mode == noscaling
                    eofs[idx] = eof(member.data[:, :, chunking[idx]]; nmodes=nmodes, center=center, align_eofs_with_mean=align_eofs_with_mean, norm_withsqrt_timedim=norm_withsqrt_timedim, weights=weights)
                else
                    eofs[idx] = scale_eof_result(eof(member.data[:, :, chunking[idx]]; nmodes=nmodes, center=center, align_eofs_with_mean=align_eofs_with_mean, norm_withsqrt_timedim=norm_withsqrt_timedim, weights=weights); scale_mode=scale_mode)
                end
            end
        else
            for idx in eachindex(chunking)
                if scale_mode == noscaling
                    eofs[idx] = eof(member.data[:, :, chunking[idx]]; nmodes=nmodes, center=center, align_eofs_with_mean=align_eofs_with_mean, norm_withsqrt_timedim=norm_withsqrt_timedim, weights=weights)
                else
                    eofs[idx] = scale_eof_result(eof(member.data[:, :, chunking[idx]]; nmodes=nmodes, center=center, align_eofs_with_mean=align_eofs_with_mean, norm_withsqrt_timedim=norm_withsqrt_timedim, weights=weights); scale_mode=scale_mode)
                end
            end
        end


        result[member.id] = eofs
        flush(stdout)
    end

    if !isnothing(saving_filepath)
        try
            save(saving_filepath, result)
        catch e
            println("Couldn't save to filepath $saving_filepath: $e")
        end

    end

    return result

end

function get_mean_of_multiple_arrays(arrays::AbstractArray...)
    sum_array = zeros(Float64, size(arrays[1]))

    for array in arrays
        sum_array += array
    end

    return sum_array ./ length(arrays)
end


function load_eof_ensemble_result(base_path, scope_id, scenario_id; sqrtscale=true, modes=5, scale_eofs::Union{Nothing,EOFScaling}=nothing, unit_scale_factor=1)::EOFEnsembleResult

    sqrt_string = sqrtscale ? "sqrtscale" : "nosqrtscale"
    ds = load(joinpath(base_path, scope_id, "eofs_$(modes)modes_$(scenario_id)_$(scope_id)_$(sqrt_string).jld2"))

    if isnothing(scale_eofs)
        ivt_piControl = convert(Vector{EOFResult}, ds["ivt_piControl"]["r1i1p1f1"])
        ps_piControl = convert(Vector{EOFResult}, ds["ps_piControl"]["r1i1p1f1"])
        ivt_eof = convert(Dict{String,Vector{EOFResult}}, ds["ivt_eof"])
        ps_eof = convert(Dict{String,Vector{EOFResult}}, ds["ps_eof"])
    else
        ivt_piControl = scale_eof_result.(convert(Vector{EOFResult}, ds["ivt_piControl"]["r1i1p1f1"]); scale_mode=scale_eofs, unit_scale_factor=unit_scale_factor)
        ps_piControl = scale_eof_result.(convert(Vector{EOFResult}, ds["ps_piControl"]["r1i1p1f1"]); scale_mode=scale_eofs, unit_scale_factor=unit_scale_factor)
        ivt_eof = Dict(member_id => scale_eof_result.(eof_results; scale_mode=scale_eofs, unit_scale_factor=unit_scale_factor) for (member_id, eof_results) in convert(Dict{String,Vector{EOFResult}}, ds["ivt_eof"]))
        ps_eof = Dict(member_id => scale_eof_result.(eof_results; scale_mode=scale_eofs, unit_scale_factor=unit_scale_factor) for (member_id, eof_results) in convert(Dict{String,Vector{EOFResult}}, ds["ps_eof"]))
    end

    return EOFEnsembleResult(
        "$scenario_id $scope_id",
        convert(Vector{UnitRange{Int}}, ds["scopes"]),
        ivt_piControl,
        ps_piControl,
        ivt_eof,
        ps_eof
    )
end

function load_eof_ensemble_common_scaling(base_path, scope_id, scenario_id, field_id::String; sqrtscale=true, modes=5, scale_eofs::Union{Nothing,EOFScaling}=nothing, align_with_first::Union{Nothing,UnitRange}=nothing)::EOFEnsemble
    (_, scaling, _, _) = get_correct_var_display(field_id)
    return load_eof_ensemble(base_path, scope_id, scenario_id, field_id; sqrtscale=sqrtscale, modes=modes, scale_eofs=scale_eofs, align_with_first=align_with_first, unit_scale_factor=scaling)
end

function load_eof_ensemble(base_path, scope_id, scenario_id, field_id::String; sqrtscale=true, modes=5, scale_eofs::Union{Nothing,EOFScaling}=nothing, align_with_first::Union{Nothing,UnitRange}=nothing, unit_scale_factor=1)::EOFEnsemble

    sqrt_string = sqrtscale ? "sqrtscale" : "nosqrtscale"
    ds = load(joinpath(base_path, field_id, scenario_id, scope_id, "$(field_id)_eofs_$(modes)modes_$(scenario_id)_$(scope_id)_$(sqrt_string).jld2"))


    function transform_eof_result(eof_res::EOFResult, results)::EOFResult

        if !isnothing(align_with_first)

            alignment_fields = []

            if align_with_first.start != 1

                for i in 1:align_with_first.start-1
                    push!(alignment_fields, nothing)
                end
            end

            for i in align_with_first

                first_field = results[1].spatial_modes[:, :, i]
                alignment_field = reshape(first_field, :)
                push!(alignment_fields, alignment_field)
            end

            eof_res = realign_modes(eof_res, alignment_fields)
        end

        if !isnothing(scale_eofs)

            eof_res = scale_eof_result(eof_res; scale_mode=scale_eofs, unit_scale_factor=unit_scale_factor)
        end

        return eof_res
    end


    piControl = convert(Vector{EOFResult}, ds["piControl"]["r1i1p1f1"])
    piControl = transform_eof_result.(piControl, Ref(piControl))
    eof_data = convert(Dict{String,Vector{EOFResult}}, ds["scenario_data"])


    foreach(eof_data) do kv
        eof_data[kv[1]] = transform_eof_result.(kv[2], Ref(eof_data[get_member_id_string(1)]))
    end

    return EOFEnsemble(
        convert(String, ds["variable_id"]),
        convert(Vector{DateTime}, ds["time"]),
        convert(Vector{Float64}, ds["lon"]),
        convert(Vector{Float64}, ds["lat"]),
        convert(Vector{UnitRange{Int}}, ds["scopes"]),
        piControl,
        eof_data
    )
end

function get_all_vertices_from_iscontours(contour_collections::Contour.ContourCollection...)::Vector{Tuple{Float64,Float64}}
    result = Tuple{Float64,Float64}[]
    for contour_collection in contour_collections
        for contour_level in levels(contour_collection)
            for contour_line in Contour.lines(contour_level)
                append!(result, [(x, y) for (x, y) in contour_line.vertices if !ismissing(x) && !ismissing(y)])
            end
        end
    end
    return result
end

function convert_curve_to_nonmissing(c::Curve2{Tuple{Union{Missing, AbstractFloat},Union{Missing, AbstractFloat}}})::Vector{Tuple{Float64,Float64}}
    
    return [(x, y) for (x, y) in c.vertices  if !ismissing(x) && !ismissing(y)]
end


function get_isocontour_lines(contour_collections::Contour.ContourCollection...)::Vector{Vector{Vector{Tuple{Float64,Float64}}}}
    res = Vector{Vector{Tuple{Float64,Float64}}}[]
    for contour_collection in contour_collections
        for contour_level in levels(contour_collection)
            push!(res, [convert_curve_to_nonmissing(c) for c in Contour.lines(contour_level)])

        end
    end

    return res
end

function get_isocontour_vertices(contour_collections::Contour.ContourCollection...)::Vector{Vector{Tuple{Float64,Float64}}}
    result = Vector{Tuple{Float64,Float64}}[]
    for contour_collection in contour_collections
        for contour_level in levels(contour_collection)
            for contour_line in Contour.lines(contour_level)
                push!(result, [(x, y) for (x, y) in contour_line.vertices if !ismissing(x) && !ismissing(y)])
            end
        end
    end
    return result
end

function get_isocontour_vertices(clines_vectors::Vector{Vector{Vector{Tuple{Float64,Float64}}}}...)::Vector{Vector{Tuple{Float64,Float64}}}
    result = Vector{Tuple{Float64,Float64}}[]

    for clines_vector in clines_vectors
        for line_vec in clines_vector
            for contour_line in line_vec
                push!(result, [(x, y) for (x, y) in contour_line])
            end
        end
    end

    return result
end

unzip(a) = map(x -> getfield.(a, x), fieldnames(eltype(a)))


function sample_along_line(line; dx=0.1)::Vector{Tuple{Float64,Float64}}
    n = length(line)
    sampled_points = []

    # Calculate total length of the line
    total_length = 0.0
    for i in 1:(n-1)
        total_length += norm([line[i+1][1] - line[i][1], line[i+1][2] - line[i][2]])
    end

    # Interpolate and sample along the line
    distance_covered = 0.0
    for i in 1:(n-1)
        (x1, y1) = line[i]
        (x2, y2) = line[i+1]
        segment_length = norm([x2 - x1, y2 - y1])

        while distance_covered < segment_length
            t = distance_covered / segment_length
            x_sample = (1 - t) * x1 + t * x2
            y_sample = (1 - t) * y1 + t * y2
            push!(sampled_points, (x_sample, y_sample))
            distance_covered += dx
        end

        # Adjust the distance covered for the next segment
        distance_covered -= segment_length
    end

    # Add the last point of the line
    push!(sampled_points, line[end])

    return sampled_points
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

function get_correct_var_display(varid::String)::Tuple{String, Float64, String, Any}

    info_dict = Dict(t[1] => t for t in [("ivt", 1.0, "kg s-1 m-1", Reverse(:vik100)), ("psl", 1/100, "hPa", :vik100), ("pr", 86400.0, "mm/month", Reverse(:vik100))])
    
    return info_dict[varid]
end
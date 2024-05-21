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

function get_member_id_string(member_nr::Int)::String
    return "r$(member_nr)i1p1f1"
end

function get_files_of_member(data_path, scenario_id, member_nr)
    return readdir(joinpath(data_path, scenario_id, get_member_id_string(member_nr)), join=true)
end

function get_data(data_path, scenario_id, member_nr; file_range_selection=:, field_id="ivt")

    file_paths = get_files_of_member(data_path, scenario_id, member_nr)

    ivt_data = Array{Float64,3}[]



    for file_path in file_paths[file_range_selection]
        ivt_chunk = ncread(file_path, field_id)
        push!(ivt_data, ivt_chunk)
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

    lons = ncread(get_files_of_member(base_path, scenarios[1], member)[1], "lon")
    lats = get_field(get_files_of_member(base_path, scenarios[1], member)[1], "lat")

    time = get_time_data(base_path, scenarios[1], member; file_range_selection=file_range_selection)

    scenarios = map(scenarios) do scenario_id

        data = get_data(base_path, scenario_id, member; file_range_selection=file_range_selection, field_id=data_field_id)

        return ScenarioData(scenario_id, data)
    end

    return TimelineData(lons, lats, time, collect(scenarios))
end

function build_ensemble_data(base_path, scenarios...; file_range_selection=:, data_field_id="ivt", member_range=1:50, silent=false, filterfun=nothing)

    lons = ncread(get_files_of_member(base_path, scenarios[1], 1)[1], "lon")
    lats = get_field(get_files_of_member(base_path, scenarios[1], 1)[1], "lat")

    result = EnsembleSimulation[]


    for scenario in scenarios
        if !silent
            println("Handling scenario $scenario ...")
        end

        time = get_time_data(base_path, scenario, 1; file_range_selection=file_range_selection)

        time_selector = isnothing(filterfun) ? Colon() : [i for i in eachindex(time) if filterfun(time[i])]
        time = time[time_selector]


        ensemble_members = map(member_range) do member_nr

            member_id = get_member_id_string(member_nr)

            data = get_data(base_path, scenario, member_nr; file_range_selection=file_range_selection, field_id=data_field_id)
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

        members[index] = EnsembleMember(pop!(member_id_set), cat([es.members[index].data for es in ensemble_simulations]..., dims=3))

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
    scale_mode=nothing,
    saving_filepath=nothing
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
        for idx in eachindex(chunking)
            eofs[idx] = eof(member.data[:, :, chunking[idx]]; nmodes=nmodes, center=center, align_eofs_with_mean=align_eofs_with_mean, norm_withsqrt_timedim=norm_withsqrt_timedim, weights=weights, scaling=scale_mode)
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


function load_eof_ensemble_result(base_path, scope_id, scenario_id; sqrtscale=false)::EOFEnsembleResult

    sqrt_string = sqrtscale ? "sqrtscale" : "nosqrtscale"
    ds = load(joinpath(base_path, scope_id, "eofs_$(scenario_id)_$(scope_id)_$(sqrt_string).jld2"))

    return EOFEnsembleResult(
        "$scenario_id $scope_id",
        convert(Vector{UnitRange{Int}}, ds["scopes"]),
        convert(Vector{EOFResult}, ds["ivt_piControl"]["r1i1p1f1"]),
        convert(Vector{EOFResult}, ds["ps_piControl"]["r1i1p1f1"]),
        convert(Dict{String,Vector{EOFResult}}, ds["ivt_eof"]),
        convert(Dict{String,Vector{EOFResult}}, ds["ps_eof"])
    )
end

function get_all_vertices_from_iscontours(contour_collections::Contour.ContourCollection...)::Vector{Tuple{Float64,Float64}}
    result = Tuple{Float64,Float64}[]
    for contour_collection in contour_collections
        for contour_level in levels(contour_struct)
            for contour_line in Contour.lines(contour_level)
                append!(result, [(x, y) for (x, y) in contour_line.vertices if !ismissing(x) && !ismissing(y)])
            end
        end
    end
    return result
end

unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))
using GLMakie
using GeoMakie
using EmpiricalOrthogonalFunctions
using NCDatasets
using NetCDF
using Dates
using BenchmarkTools
using Statistics



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

struct EOFResult
    spatial_modes::Array{Float64, 3}
    temporal_modes::Array{Float64, 2}
    modes_variability::Array{Float64, 1}
end    

function get_files_of_member(data_path, scenario_id, member_nr)
    return readdir(joinpath(data_path, scenario_id, "r$(member_nr)i1p1f1"), join=true)
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

function filter_by_date(fun, timeline_data::TimelineData)::TimelineData

    time_indices = [i for i in eachindex(timeline_data.time) if fun(timeline_data.time[i])]

    transformed_scenarios = [ScenarioData(ds.name, ds.data[:, :, time_indices]) for ds in timeline_data.datasets]

    return TimelineData(timeline_data.lons, timeline_data.lats, timeline_data.time[time_indices], transformed_scenarios)
end

function align_with_field!(field, alignment_field; dims = 1)

    for slice in eachslice(field, dims = dims)
        scalar_product = sum(slice .* alignment_field)
        if scalar_product < 0
            slice *= -1
        end
    end
end

function align_with_field(field, alignment_field; dims = 1)

    function handle_slice(slice)
        scalar_product = sum(slice .* alignment_field)
        if scalar_product < 0
            return slice * -1
        else
            return slice
        end
    end

    return cat([handle_slice(slice) for slice in eachslice(field, dims = dims)]..., dims = 2)
end

function get_eof_of_datachunk(datachunk; nmodes = nothing, center = true, alignment_field = nothing)::EOFResult

    eof = EmpiricalOrthogonalFunction(datachunk; timedim=3, center = center)

    if isnothing(nmodes)
        nmodes = size(datachunk, 3)
    end

    temporalsignal = pcs(eof)
    spatialsignal = eofs(eof)

    if !isnothing(alignment_field)
        spatialsignal = align_with_field(spatialsignal, alignment_field; dims = 1)
    end

    modes_variability = eof.eigenvals ./ sum(eof.eigenvals) * 100
    return EOFResult(reshape(spatialsignal, (size(datachunk)[1:2]..., :))[:, :, 1:nmodes], temporalsignal[:, 1:nmodes], modes_variability[1:nmodes])
end



function get_spatial_modes(datachunk; nmodes = nothing, center = true)
    eof = EmpiricalOrthogonalFunction(datachunk; timedim=3, center = center)

    if isnothing(nmodes)
        nmodes = size(datachunk, 3)
    end

    temporalsignal = pcs(eof)
    spatialsignal = reshape(eofs(eof), (size(datachunk)[1:2]..., :))
    modes_variability = eof.eigenvals ./ sum(eof.eigenvals) * 100
    return spatialsignal[:, :, 1:nmodes], temporalsignal[:, 1:nmodes], modes_variability[1:nmodes]
end


function get_sliding_time_scopes_by_threshold(time_data, n, threshold = Day(100))
    
    all_starting_points = [1]

    for i in 1:length(time_data)-1 
        if time_data[i+1] - time_data[i] > threshold
            push!(all_starting_points, i+1)
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

function centralize_tl_data(tl_data; dims = 3)

    normalized_ds = ScenarioData[]
    for ds in tl_data.datasets
        res = ds.data .- mean(ds.data, dims = dims)
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
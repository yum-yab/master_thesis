using GLMakie
using GeoMakie
using EmpiricalOrthogonalFunctions
using NCDatasets
using Dates
using BenchmarkTools
using Statistics
using StatsBase
using Colors
using Contour
using NumericalIntegration
using ProgressMeter

include("utils.jl")
include("ensemble_hexbin.jl")


function build_figure(
    data,
    lon_bounds::Tuple{<:Number,<:Number},
    lat_bounds::Tuple{<:Number,<:Number};
    title::String="",
    colormap=:viridis,
    shading=Makie.automatic,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing
)::Makie.Figure

    fig = isnothing(resolution) ? Figure(fontsize=12) : Figure(size=resolution, fontsize=12)

    minval, maxval = extrema(data)

    ga = local_geoaxis_creation!(fig, lon_bounds, lat_bounds; title=title)
    surface!(ga, lon_bounds[1] .. lon_bounds[2], lat_bounds[1] .. lat_bounds[2], data; shading=shading, colormap=colormap)
    lines!(ga, GeoMakie.coastlines(); color=:white, transformation=(; translation=(0, 0, 1000)))
    Colorbar(fig[1, 2], limits=(minval, maxval), colormap=colormap)
    return fig
end

function animate_timeline(
    data::TimelineData,
    filename::String;
    framerate::Int=30,
    colormap=:viridis,
    coastline_color=:white,
    shading=Makie.automatic,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing
)::Makie.Figure

    fig = isnothing(resolution) ? Figure(fontsize=12) : Figure(size=resolution, fontsize=12)

    time_index = Observable(1)

    axis = Dict()

    lon_bounds = extrema(data.lons)
    lat_bounds = extrema(data.lats)

    all_extrema = extrema.([d.data for d in data.datasets])

    min_val, max_val = reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), all_extrema)
    Colorbar(fig[2, 1:2], limits=(min_val, max_val), colormap=colormap, vertical=false, label="IVT (norm) in kg s^-1 m^-1")

    for (i, dataset) in enumerate(data.datasets)
        axis[dataset.name] = local_geoaxis_creation!(fig, lon_bounds, lat_bounds; title=@lift("$(dataset.name) $(Dates.format(data.time[$time_index], "mm/yyyy"))"), figure_col=i)
    end

    for dataset in data.datasets

        array_slice = @lift(dataset.data[:, :, $time_index])

        surface!(axis[dataset.name], lon_bounds[1] .. lon_bounds[2], lat_bounds[1] .. lat_bounds[2], array_slice; shading=shading, colormap=colormap, colorrange=(min_val, max_val))
        lines!(axis[dataset.name], GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))
    end

    timestamps = eachindex(data.time)
    record(fig, filename, timestamps; framerate=framerate) do t
        time_index[] = t
    end


    return fig
end

function get_extrema_of_arrays(datasets::AbstractArray{<:Number}...)

    all_extrema = Tuple{Float64,Float64}[]

    for ds in datasets
        push!(all_extrema, extrema(ds))
    end

    return reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), all_extrema)

end

function show_eof_modes_of_timeline(
    data::TimelineData,
    eof_data::Dict{String,Vector{EOFResult}},
    all_scopes,
    filename::String;
    framerate::Int=30,
    colormap=:viridis,
    coastline_color=:white,
    shading=Makie.automatic,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing,
    fontsize::Int=12,
    whole_figure_title=nothing,
)::Makie.Figure

    nmodes = length(eof_data[data.datasets[1].name][1].modes_variability)
    fig = isnothing(resolution) ? Figure(fontsize=fontsize) : Figure(size=resolution, fontsize=fontsize)

    if !isnothing(whole_figure_title)
        Label(fig[0, 1:nmodes], whole_figure_title, fontsize=round(1.7 * fontsize))
    end

    current_scope_index = Observable(1)

    axis = Dict()

    lon_bounds = extrema(data.lons)
    lat_bounds = extrema(data.lats)

    eof_extremas = Tuple{Float64,Float64}[]


    for (i, dataset) in enumerate(data.datasets)
        all_modes_axis = []


        axis[dataset.name] = all_modes_axis
        append!(eof_extremas, [extrema(res.spatial_modes) for res in eof_data[dataset.name]])


        for j in 1:nmodes
            title = @lift("$(dataset.name) $(year(data.time[all_scopes[$current_scope_index].start]))-$(year(data.time[all_scopes[$current_scope_index].stop])) Mode $j Variability $(round(eof_data[dataset.name][$current_scope_index].modes_variability[j], digits = 2))")
            push!(all_modes_axis, local_geoaxis_creation!(fig, lon_bounds, lat_bounds; title=title, figure_row=i, figure_col=j))
        end
    end

    min_val, max_val = reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), eof_extremas)

    Colorbar(fig[length(data.datasets)+1, 1:nmodes], limits=(min_val, max_val), colormap=colormap, vertical=false, label="EOF values")

    for (i, dataset) in enumerate(data.datasets)

        for mode in 1:nmodes
            surface!(axis[dataset.name][mode], lon_bounds[1] .. lon_bounds[2], lat_bounds[1] .. lat_bounds[2], @lift(eof_data[dataset.name][$current_scope_index].spatial_modes[:, :, mode]); shading=shading, colormap=colormap, colorrange=(min_val, max_val))
            lines!(axis[dataset.name][mode], GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))
        end
    end

    record(fig, filename, eachindex(all_scopes); framerate=framerate) do t
        current_scope_index[] = t
    end


    return fig
end

function compare_truth_with_tldata(
    truth::TimelineData,
    data::TimelineData,
    filename::String;
    framerate::Int=30,
    colormap=:viridis,
    coastline_color=:white,
    shading=Makie.automatic,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing
)::Makie.Figure

    fig = isnothing(resolution) ? Figure(fontsize=12) : Figure(size=resolution, fontsize=12)

    time_index = Observable(1)


    axis = Dict()

    lon_bounds = extrema(data.lons)
    lat_bounds = extrema(data.lats)

    all_extrema = extrema.([[d.data for d in data.datasets]; [d.data for d in truth.datasets]])

    min_val, max_val = reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), all_extrema)

    axis["truth"] = local_geoaxis_creation!(fig, lon_bounds, lat_bounds; title=@lift("Reference ssp126 $(Dates.format(truth.time[$time_index], "mm/yyyy"))"))

    for (i, dataset) in enumerate(data.datasets)
        axis[dataset.name] = local_geoaxis_creation!(fig, lon_bounds, lat_bounds; title=@lift("$(dataset.name) $(Dates.format(data.time[$time_index], "mm/yyyy"))"), figure_col=i + 1)
    end

    Colorbar(fig[2, 1:length(axis)], limits=(min_val, max_val), colormap=colormap, vertical=false, label="IVT (norm) in kg s^-1 m^-1")

    truth_slice = @lift(truth.datasets[1].data[:, :, $time_index])

    surface!(axis["truth"], lon_bounds[1] .. lon_bounds[2], lat_bounds[1] .. lat_bounds[2], truth_slice; shading=shading, colormap=colormap, colorrange=(min_val, max_val))
    lines!(axis["truth"], GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

    for dataset in data.datasets

        array_slice = @lift(dataset.data[:, :, $time_index])

        surface!(axis[dataset.name], lon_bounds[1] .. lon_bounds[2], lat_bounds[1] .. lat_bounds[2], array_slice; shading=shading, colormap=colormap, colorrange=(min_val, max_val))
        lines!(axis[dataset.name], GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))
    end

    timestamps = eachindex(truth.time)
    record(fig, filename, timestamps; framerate=framerate) do t
        time_index[] = t
    end


    return fig
end

function show_eof_and_pc_modes_of_timeline(
    data::TimelineData,
    eof_data::Dict{String,Vector{EOFResult}},
    all_scopes,
    filename::String;
    framerate::Int=30,
    colormap=:viridis,
    coastline_color=:white,
    shading=Makie.automatic,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing,
    fontsize::Int=12,
    whole_figure_title=nothing,
    surface_alpha_value=1.0,
    additional_eof_data::Dict{String,Vector{EOFResult}}=nothing,
    additional_data_label="Additional Data Label",
)::Makie.Figure

    nmodes = length(eof_data[data.datasets[1].name][1].modes_variability)

    additional_data_modes = isnothing(additional_eof_data) ? 0 : length(additional_eof_data[collect(keys(additional_eof_data))[1]][1].modes_variability)



    fig = isnothing(resolution) ? Figure(fontsize=fontsize) : Figure(size=resolution, fontsize=fontsize)

    if !isnothing(whole_figure_title)
        Label(fig[0, 1:nmodes+1], whole_figure_title, fontsize=round(1.7 * fontsize))
    end

    current_scope_index = Observable(1)

    eof_data_axis = Dict()

    pc_axis = Dict()

    mean_axis = Dict()

    lonmin, lonmax = extrema(data.lons)
    latmin, latmax = extrema(data.lats)

    eof_extremas = Tuple{Float64,Float64}[]

    pcs_extremas = Tuple{Float64,Float64}[]

    for (i, dataset) in enumerate(data.datasets)
        all_modes_axis = []


        eof_data_axis[dataset.name] = all_modes_axis
        pc_axis[dataset.name] = Axis(fig[2*i, 1:nmodes+1], limits=((nothing, nothing), (-1, 1)))
        for res in eof_data[dataset.name]
            push!(eof_extremas, extrema(res.spatial_modes))
            push!(pcs_extremas, extrema(res.temporal_modes))
        end


        for j in 1:nmodes
            title = @lift("$(dataset.name) $(year(data.time[all_scopes[$current_scope_index].start]))-$(year(data.time[all_scopes[$current_scope_index].stop])) Mode $j Variability $(round(eof_data[dataset.name][$current_scope_index].modes_variability[j], digits = 2)) %")
            push!(all_modes_axis, local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=title, figure_row=2 * (i - 1) + 1, figure_col=j))
        end
        # also create axis for mean map
        mean_axis[dataset.name] = local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=@lift("$(dataset.name) $(year(data.time[all_scopes[$current_scope_index].start]))-$(year(data.time[all_scopes[$current_scope_index].stop])) Mean map"), figure_row=2 * (i - 1) + 1, figure_col=nmodes + 1)
    end

    spatial_min_val, spatial_max_val = reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), eof_extremas)

    println("Data extremas: $spatial_min_val $spatial_max_val")

    Colorbar(fig[2*length(data.datasets)+1, 1:nmodes+1], limits=(spatial_min_val, spatial_max_val), colormap=colormap, vertical=false, label="EOF values")

    for (i, dataset) in enumerate(data.datasets)

        dataset_mean = @lift(dropdims(mean(dataset.data[:, :, all_scopes[$current_scope_index]], dims=3), dims=3))

        surface!(mean_axis[dataset.name], lonmin .. lonmax, latmin .. latmax, dataset_mean; shading=shading, colormap=colormap, colorrange=(spatial_min_val, spatial_max_val), alpha=surface_alpha_value, overdraw=false, transformation=(; translation=(0, 0, -1000)))
        lines!(mean_axis[dataset.name], GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))
        for mode in 1:nmodes
            surface!(eof_data_axis[dataset.name][mode], lonmin .. lonmax, latmin .. latmax, @lift(eof_data[dataset.name][$current_scope_index].spatial_modes[:, :, mode]); shading=shading, colormap=colormap, colorrange=(spatial_min_val, spatial_max_val), alpha=surface_alpha_value, overdraw=false, transformation=(; translation=(0, 0, -1000)))
            lines!(eof_data_axis[dataset.name][mode], GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

            #positions = @lift([(i, eof_data[dataset.name][$current_scope_index].temporal_modes[i, mode]) for i in all_scopes[$current_scope_index]])
            positions = @lift(collect(zip(eachindex(all_scopes[$current_scope_index]), standardize(UnitRangeTransform, eof_data[dataset.name][$current_scope_index].temporal_modes[:, mode]; unit=false))))

            #lines!(pc_axis[dataset.name], 1..30, @lift(eof_data[dataset.name][$current_scope_index].temporal_modes[:, mode]); label = "Mode $mode")

            # draw pcs lines 
            lines!(pc_axis[dataset.name], positions; label="Mode $mode")

            #draw zero contour plot
            contour!(eof_data_axis[dataset.name][mode], lonmin .. lonmax, latmin .. latmax, @lift(eof_data[dataset.name][$current_scope_index].spatial_modes[:, :, mode]), levels=[0], color=:red, transformation=(; translation=(0, 0, 1100)))

            # draw the additional data
            if !isnothing(additional_eof_data) && mode in 1:additional_data_modes && dataset.name in keys(additional_eof_data)

                positions_additional_data = @lift(collect(zip(eachindex(all_scopes[$current_scope_index]), standardize(UnitRangeTransform, additional_eof_data[dataset.name][$current_scope_index].temporal_modes[:, mode]; unit=false))))
                lines!(pc_axis[dataset.name], positions_additional_data; label="$additional_data_label Mode $mode")
            end

        end
        axislegend(pc_axis[dataset.name])
    end

    # time_data = string.(data.time)
    # dates_slice = @lift(range(all_scopes[$current_scope_index].start, all_scopes[$current_scope_index].stop, step=4))
    vals_per_season = all_scopes[3].start - all_scopes[2].start
    vals_per_scope = length(all_scopes[1])
    seasonslice = range(1, vals_per_scope, step=vals_per_season)

    for (_, axis) in pc_axis
        axis.xticks = (seasonslice, ["Winter $i" for i in eachindex(seasonslice)])
        axis.xticklabelrotation = π / 4
        axis.xticklabelalign = (:right, :center)
    end

    record(fig, filename, eachindex(all_scopes); framerate=framerate) do t
        current_scope_index[] = t
    end


    return fig
end

function show_eof_and_pc_modes_of_member(
    data::EnsembleSimulation,
    eof_data::Dict{String,Vector{EOFResult}},
    all_scopes,
    member,
    filename::String;
    framerate::Int=30,
    colormap=:viridis,
    coastline_color=:white,
    shading=Makie.automatic,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing,
    fontsize::Int=12,
    whole_figure_title=nothing,
    surface_alpha_value=1.0,
    draw_contours=true,
    nmodes=nothing,
)::Makie.Figure


    member_identifier = get_member_id_string(member)
    nmodes = isnothing(nmodes) ? length(eof_data[member_identifier][1].singularvals) : nmodes


    fig = isnothing(resolution) ? Figure(fontsize=fontsize) : Figure(size=resolution, fontsize=fontsize)

    current_scope_index = Observable(1)


    if !isnothing(whole_figure_title)
        Label(fig[0, 1:nmodes+1], @lift("$(whole_figure_title) $(member_identifier) $(year(data.time[all_scopes[$current_scope_index].start]))-$(year(data.time[all_scopes[$current_scope_index].stop]))"), fontsize=round(1.7 * fontsize))
    end



    pc_axis = Axis(fig[3, 1:nmodes+1], limits=((nothing, nothing), (-1, 1)))

    lonmin, lonmax = extrema(data.lons)
    latmin, latmax = extrema(data.lats)


    spatial_min_val, spatial_max_val = reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), [extrema(res.spatial_modes) for res in eof_data[member_identifier]])

    data_min_val, data_max_val = reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), [extrema(mean(data.members[member].data[:, :, scope], dims=3)) for scope in all_scopes])
    # spatial_min_val = -1 * spatial_max_val
    println("temporal extremata: $(extrema(reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), [extrema(res.temporal_modes[:, 1:nmodes]) for res in eof_data[member_identifier]])))")
    println("EOF extremas: $spatial_min_val $spatial_max_val")
    println("Data mean extremas: $data_min_val $data_max_val")

    upper_limit = max(data_max_val, spatial_max_val)

    lower_limit = -1 * upper_limit


    mode_axes = []

    for modenr in 1:nmodes
        title = @lift("Mode $modenr Variability $(round(get_modes_variability(eof_data[member_identifier][$current_scope_index])[modenr] * 100, digits = 2)) %")
        push!(mode_axes, local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=title, figure_row=1, figure_col=modenr))
    end

    mean_axis = local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=@lift("$(data.id) $(year(data.time[all_scopes[$current_scope_index].start]))-$(year(data.time[all_scopes[$current_scope_index].stop])) Mean map"), figure_row=1, figure_col=nmodes + 1)

    dataset_mean_observable = @lift(dropdims(mean(data.members[member].data[:, :, all_scopes[$current_scope_index]], dims=3), dims=3))

    surface!(mean_axis, lonmin .. lonmax, latmin .. latmax, dataset_mean_observable; shading=shading, colormap=colormap, colorrange=(lower_limit, upper_limit), alpha=surface_alpha_value, overdraw=false, transformation=(; translation=(0, 0, -1000)))
    lines!(mean_axis, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))


    for mode in 1:nmodes
        surface!(mode_axes[mode], lonmin .. lonmax, latmin .. latmax, @lift(eof_data[member_identifier][$current_scope_index].spatial_modes[:, :, mode]), shading=shading, colormap=colormap, colorrange=(lower_limit, upper_limit), alpha=surface_alpha_value, overdraw=false, transformation=(; translation=(0, 0, -1 * spatial_max_val)))
        lines!(mode_axes[mode], GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))
        positions = @lift(collect(zip(eachindex(all_scopes[$current_scope_index]), standardize(UnitRangeTransform, eof_data[member_identifier][$current_scope_index].temporal_modes[:, mode]; unit=false))))
        lines!(pc_axis, positions; label="Mode $mode")
    end
    axislegend(pc_axis)

    if draw_contours

        for mode in 1:nmodes
            contour!(mode_axes[mode], lonmin .. lonmax, latmin .. latmax, @lift(eof_data[member_identifier][$current_scope_index].spatial_modes[:, :, mode]), levels=[0], color=:red, transformation=(; translation=(0, 0, 1100)))
        end
    end



    Colorbar(fig[2, 1:nmodes+1], limits=(lower_limit, upper_limit), colormap=colormap, vertical=false, label="EOF values")

    # time_data = string.(data.time)
    # dates_slice = @lift(range(all_scopes[$current_scope_index].start, all_scopes[$current_scope_index].stop, step=4))
    vals_per_season = all_scopes[3].start - all_scopes[2].start
    vals_per_scope = length(all_scopes[1])
    seasonslice = range(1, vals_per_scope, step=vals_per_season)

    pc_axis.xticks = (seasonslice, ["Winter $i" for i in eachindex(seasonslice)])
    pc_axis.xticklabelrotation = π / 4
    pc_axis.xticklabelalign = (:right, :center)


    record(fig, filename, eachindex(all_scopes); framerate=framerate) do t
        current_scope_index[] = t
    end


    return fig
end



function show_eof_and_pc_modes_of_ensemble(
    data::EnsembleSimulation,
    eof_data::Dict{String,Vector{EOFResult}},
    all_scopes,
    filename::String;
    framerate::Int=30,
    colormap=:viridis,
    coastline_color=:white,
    shading=Makie.automatic,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing,
    fontsize::Int=12,
    whole_figure_title=nothing,
    surface_alpha_value=1.0,
    draw_contours=true
)::Makie.Figure

    nmodes = length(eof_data[get_member_id_string(1)][1].singularvals)


    fig = isnothing(resolution) ? Figure(fontsize=fontsize) : Figure(size=resolution, fontsize=fontsize)

    if !isnothing(whole_figure_title)
        Label(fig[0, 1:nmodes+1], whole_figure_title, fontsize=round(1.7 * fontsize))
    end

    current_scope_index = Observable(1)


    pc_axis = Axis(fig[3, 1:nmodes+1], limits=((nothing, nothing), (-1, 1)))

    lonmin, lonmax = extrema(data.lons)
    latmin, latmax = extrema(data.lats)

    eof_extremas = Tuple{Float64,Float64}[]



    for (_, eofres_array) in eof_data
        append!(eof_extremas, [extrema(res.spatial_modes) for res in eofres_array])
    end


    spatial_min_val, spatial_max_val = reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), eof_extremas)

    data_min_val, data_max_val = extrema(get_mean_of_multiple_arrays([member.data for member in data.members]...))
    # spatial_min_val = -1 * spatial_max_val
    println("EOF extremas: $spatial_min_val $spatial_max_val")
    println("Data extremas: $data_min_val $data_max_val")

    upper_limit = max(data_max_val, spatial_max_val)

    lower_limit = -1 * upper_limit


    mode_axes = []

    for modenr in 1:nmodes
        title = @lift("$(data.id) $(year(data.time[all_scopes[$current_scope_index].start]))-$(year(data.time[all_scopes[$current_scope_index].stop])) Mode $modenr Variability $(round(mean([get_modes_variability(resarray[$current_scope_index])[modenr] for (_, resarray) in eof_data]), digits = 2)) %")
        push!(mode_axes, local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=title, figure_row=1, figure_col=modenr))
    end

    mean_axis = local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=@lift("$(data.id) $(year(data.time[all_scopes[$current_scope_index].start]))-$(year(data.time[all_scopes[$current_scope_index].stop])) Mean map"), figure_row=1, figure_col=nmodes + 1)

    dataset_mean_observable = @lift(dropdims(mean(get_mean_of_multiple_arrays([member.data[:, :, all_scopes[$current_scope_index]] for member in data.members]...), dims=3), dims=3))

    surface!(mean_axis, lonmin .. lonmax, latmin .. latmax, dataset_mean_observable; shading=shading, colormap=colormap, colorrange=(lower_limit, upper_limit), alpha=surface_alpha_value, overdraw=false, transformation=(; translation=(0, 0, -1000)))
    lines!(mean_axis, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))


    for mode in 1:nmodes
        surface!(mode_axes[mode], lonmin .. lonmax, latmin .. latmax, @lift(get_mean_of_multiple_arrays([res_array[$current_scope_index].spatial_modes[:, :, mode] for (_, res_array) in eof_data]...)); shading=shading, colormap=colormap, colorrange=(lower_limit, upper_limit), alpha=surface_alpha_value, overdraw=false, transformation=(; translation=(0, 0, -1 * spatial_max_val)))
        lines!(mode_axes[mode], GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))
        positions = @lift(collect(zip(eachindex(all_scopes[$current_scope_index]), standardize(UnitRangeTransform, get_mean_of_multiple_arrays([res_array[$current_scope_index].temporal_modes[:, mode] for (_, res_array) in eof_data]...); unit=false))))
        lines!(pc_axis, positions; label="Mode $mode")
    end
    axislegend(pc_axis)

    if draw_contours
        contour_colors = distinguishable_colors(50)

        for (i, (_, eofres_array)) in enumerate(eof_data)
            for mode in 1:nmodes
                contour!(mode_axes[mode], lonmin .. lonmax, latmin .. latmax, @lift(eofres_array[$current_scope_index].spatial_modes[:, :, mode]), levels=[0], color=contour_colors[i], transformation=(; translation=(0, 0, 1100)))
            end
        end
    end



    Colorbar(fig[2, 1:nmodes+1], limits=(lower_limit, upper_limit), colormap=colormap, vertical=false, label="EOF values")

    # time_data = string.(data.time)
    # dates_slice = @lift(range(all_scopes[$current_scope_index].start, all_scopes[$current_scope_index].stop, step=4))
    vals_per_season = all_scopes[3].start - all_scopes[2].start
    vals_per_scope = length(all_scopes[1])
    seasonslice = range(1, vals_per_scope, step=vals_per_season)

    pc_axis.xticks = (seasonslice, ["Winter $i" for i in eachindex(seasonslice)])
    pc_axis.xticklabelrotation = π / 4
    pc_axis.xticklabelalign = (:right, :center)


    record(fig, filename, eachindex(all_scopes); framerate=framerate) do t
        current_scope_index[] = t
    end


    return fig
end


function compare_data_reconstruction(
    ensemble_simulation::EnsembleSimulation,
    eof_result::EOFEnsembleResult,
    member::Int,
    filename::String;
    scope_index=1,
    field=:ivt,
    framerate::Int=30,
    data_colormap=:vik100,
    coastline_color=:white,
    shading=Makie.automatic,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing,
    fontsize::Int=12,
    whole_figure_title=nothing,
    filter_by_temporal_mode_fun=nothing
)

    ensemble_eofs, _ = pick_correct_field(eof_result, field)

    current_index = Observable(1)
    lonmin, lonmax = extrema(ensemble_simulation.lons)
    latmin, latmax = extrema(ensemble_simulation.lats)

    my_scope = eof_result.scopes[scope_index]

    scope_timeaxis = ensemble_simulation.time[my_scope]

    fig = isnothing(resolution) ? Figure(fontsize=fontsize) : Figure(size=resolution, fontsize=fontsize)


    simulation_data = ensemble_simulation.members[member].data[:, :, my_scope]

    eof_res = scale_eof_result(ensemble_eofs[get_member_id_string(member)][scope_index])

    data_limits = extrema(simulation_data)


    limits = (-1 * data_limits[2], data_limits[2])

    _, eof_maxval = extrema(eof_res.spatial_modes)

    eof_limits = (-1 * eof_maxval, eof_maxval)

    nmodes = length(eof_res.singularvals)

    chosen_indices = isnothing(filter_by_temporal_mode_fun) ? eachindex(my_scope) : filter_by_temporal_mode_fun(eof_res.temporal_modes)

    if !isnothing(whole_figure_title)
        Label(fig[0, 1:nmodes+1], whole_figure_title, fontsize=round(1.5 * fontsize))
    end

    month_year_format = DateFormat("yyyy-mm")


    orig_data_ax = local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=@lift("$(ensemble_simulation.id) $(Dates.format(scope_timeaxis[$current_index], month_year_format)) Original Data"), figure_row=3, figure_col=1)

    surface!(
        orig_data_ax,
        lonmin .. lonmax,
        latmin .. latmax,
        @lift(simulation_data[:, :, $current_index]);
        shading=shading,
        colormap=data_colormap,
        colorrange=limits,
        overdraw=false,
        transformation=(; translation=(0, 0, -1 * data_limits[2]))
    )
    lines!(orig_data_ax, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

    influence_axis = Axis(fig[1, 1], limits=((nothing, nothing), (-1, 1)), title="Coefficient of Modes (standardized)")

    standardized_temp_modes = cat([standardize(UnitRangeTransform, eof_res.temporal_modes[:, i]; unit=false) for i in 1:nmodes]..., dims=2)

    positions_observable = @lift([(i, standardized_temp_modes[$current_index, i]) for i in 1:nmodes])

    barplot!(influence_axis, positions_observable, color=@lift(unzip($positions_observable)[1]))


    for modenum in 1:nmodes

        truncated_eof_res = truncate_eof_result(eof_res, modenum)
        var_encoded = sum(get_modes_variability(truncated_eof_res) * 100)

        modes_var = get_modes_variability(eof_res) * 100

        reconstructed_eof_res = reconstruct_data(truncated_eof_res;)

        println("Min and max error of reconstruction of Mode $modenum: $(extrema(abs.(simulation_data .- reconstructed_eof_res)))")

        mode_ax = local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title="Mode $modenum Var: $(modes_var[modenum]) %", figure_row=1, figure_col=1 + modenum)
        reconstructed_ax = local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title="Reconstruction Modes 1:$(modenum) Var: $(var_encoded) %", figure_row=3, figure_col=1 + modenum)

        surface!(
            mode_ax,
            lonmin .. lonmax,
            latmin .. latmax,
            spatial_modes[:, :, modenum];
            shading=shading,
            colormap=data_colormap,
            colorrange=limits,
            overdraw=false,
            transformation=(; translation=(0, 0, -1 * eof_limits[2]))
        )
        lines!(mode_ax, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

        surface!(
            reconstructed_ax,
            lonmin .. lonmax,
            latmin .. latmax,
            @lift(reconstructed_eof_res[:, :, $current_index]);
            shading=shading,
            colormap=data_colormap,
            colorrange=limits,
            overdraw=false,
            transformation=(; translation=(0, 0, -1 * data_limits[2]))
        )
        lines!(reconstructed_ax, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))


    end

    Colorbar(fig[2, 1:6], limits=limits, colormap=data_colormap, vertical=false, label="IVT Values")
    # Colorbar(fig[2, 5:6], limits=eof_limits, colormap=eof_colormap, vertical=false, label="IVT anomalies")



    record(fig, filename, chosen_indices; framerate=framerate) do t
        current_index[] = t
    end


    return fig





end


function draw_reconstruction_field!(
    figure,
    row_id,
    ensemble_simulation::EnsembleSimulation,
    eof_result::EOFEnsemble,
    member::Int,
    field::Symbol,
    modes::Int,
    current_index;
    data_label::String="",
    scope_index=1,
    data_colormap=:vik100,
    coastline_color=:white,
    shading=NoShading,
    fontsize::Int=12,
    whole_figure_title=nothing,
    display_mode=:anomaly,
    factor_for_data_transform=1
)



    lonmin, lonmax = extrema(ensemble_simulation.lons)
    latmin, latmax = extrema(ensemble_simulation.lats)

    my_scope = eof_result.scopes[scope_index]

    scope_timeaxis = ensemble_simulation.time[my_scope]

    if display_mode == :anomaly
        simulation_data = (ensemble_simulation.members[member].data[:, :, my_scope] .- mean(ensemble_simulation.members[member].data[:, :, my_scope], dims=3)) .* factor_for_data_transform
        old_mean = nothing
        data_name = "Anomaly Map"
    elseif display_mode == :reconstruction
        simulation_data = ensemble_simulation.members[member].data[:, :, my_scope]
        old_mean = mean(ensemble_simulation.members[member].data[:, :, my_scope], dims=3)
        data_name = "Original Data"
    else
        ArgumentError("Only :anomaly and :reconstruction are supported")
    end


    eof_res = eof_result.ensemble[get_member_id_string(member)][scope_index]

    spatial_modes = eof_res.spatial_modes .* factor_for_data_transform

    data_min, data_max = extrema(simulation_data)

    eof_minval, eof_maxval = extrema(spatial_modes[:, :, 1:modes])

    upper_limit = max(eof_maxval, data_max)

    limits = (-1 * eof_maxval * 3, eof_maxval * 3)

    day_format = DateFormat("yyyy-mm-dd")

    println("extrema of $data_label EOF: $(extrema(spatial_modes))")

    if !isnothing(whole_figure_title)
        Label(figure[row_id, 1:modes+1], @lift("$(whole_figure_title) $(Dates.format(scope_timeaxis[$current_index], day_format))"), fontsize=round(1.7 * fontsize))
    end

    month_year_format = DateFormat("yyyy-mm")


    orig_data_ax = local_geoaxis_creation!(figure, (lonmin, lonmax), (latmin, latmax); title=@lift("$(ensemble_simulation.id) $(Dates.format(scope_timeaxis[$current_index], month_year_format)) $data_name"), figure_row=row_id + 3, figure_col=1)

    surface!(
        orig_data_ax,
        lonmin .. lonmax,
        latmin .. latmax,
        @lift(simulation_data[:, :, $current_index]);
        shading=shading,
        colormap=data_colormap,
        colorrange=limits,
        overdraw=false,
        transformation=(; translation=(0, 0, -1 * limits[2]))
    )
    lines!(orig_data_ax, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

    influence_axis = Axis(figure[row_id+1, 1], limits=((nothing, nothing), (-1, 1)), title="Coefficient of Modes (standardized)")

    standardized_temp_modes = cat([standardize(UnitRangeTransform, eof_res.temporal_modes[:, i]; unit=false) for i in 1:modes]..., dims=2)

    positions_observable = @lift([(i, standardized_temp_modes[$current_index, i]) for i in 1:modes])

    barplot!(influence_axis, positions_observable, color=@lift(unzip($positions_observable)[1]))


    for modenum in 1:modes

        truncated_eof_res = truncate_eof_result(eof_res, modenum)
        var_encoded = sum(get_modes_variability(truncated_eof_res) * 100)

        modes_var = get_modes_variability(eof_res) * 100

        reconstructed_eof_res = reconstruct_data(truncated_eof_res; original_data_timemean=old_mean) .* factor_for_data_transform

        println("Mean Squared error of reconstruction of Mode $modenum: $(mean((simulation_data .- reconstructed_eof_res).^2))")

        mode_ax = local_geoaxis_creation!(figure, (lonmin, lonmax), (latmin, latmax); title="Mode $modenum Var: $(modes_var[modenum]) %", figure_row=row_id + 1, figure_col=1 + modenum)
        reconstructed_ax = local_geoaxis_creation!(figure, (lonmin, lonmax), (latmin, latmax); title="Reconstruction Modes 1:$(modenum) Var: $(var_encoded) %", figure_row=row_id + 3, figure_col=1 + modenum)

        surface!(
            mode_ax,
            lonmin .. lonmax,
            latmin .. latmax,
            spatial_modes[:, :, modenum];
            shading=shading,
            colormap=data_colormap,
            colorrange=limits,
            overdraw=false,
            transformation=(; translation=(0, 0, -1 * limits[2]))
        )
        lines!(mode_ax, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

        surface!(
            reconstructed_ax,
            lonmin .. lonmax,
            latmin .. latmax,
            @lift(reconstructed_eof_res[:, :, $current_index]);
            shading=shading,
            colormap=data_colormap,
            colorrange=limits,
            overdraw=false,
            transformation=(; translation=(0, 0, -1 * limits[2]))
        )
        lines!(reconstructed_ax, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))


    end

    Colorbar(figure[row_id+2, 1:1+modes], limits=limits, colormap=data_colormap, vertical=false, label=data_label)
    # Colorbar(fig[2, 5:6], limits=eof_limits, colormap=eof_colormap, vertical=false, label="IVT anomalies")
end

function compare_ivt_ps_reconstruction(
    ivt_ensemble_simulation::EnsembleSimulation,
    ps_ensemble_simulation::EnsembleSimulation,
    eof_result::EOFEnsembleResult,
    member::Int,
    modes::Int,
    filename::String;
    scope_index=1,
    framerate::Int=30,
    data_colormap=:vik100,
    coastline_color=:white,
    shading=NoShading,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing,
    fontsize::Int=12,
    whole_figure_title=nothing,
    filter_by_temporal_mode_fun=nothing
)


    current_index = Observable(1)

    my_scope = eof_result.scopes[scope_index]

    fig = isnothing(resolution) ? Figure(fontsize=fontsize) : Figure(size=resolution, fontsize=fontsize)


    chosen_indices = isnothing(filter_by_temporal_mode_fun) ? eachindex(my_scope) : filter_by_temporal_mode_fun(eof_result.ivt_ensemble_eofs[get_member_id_string(member)][scope_index])

    if !isnothing(whole_figure_title)
        Label(fig[0, 1:modes+1], whole_figure_title, fontsize=round(1.7 * fontsize))
    end


    draw_reconstruction_field!(
        fig,
        1,
        ivt_ensemble_simulation,
        eof_result,
        member,
        :ivt,
        modes,
        current_index;
        data_label="IVT Values",
        scope_index=scope_index,
        data_colormap=data_colormap,
        coastline_color=coastline_color,
        shading=shading,
        fontsize=fontsize,
        whole_figure_title="Anomalies of IVT fields",
        display_mode=:anomaly
    )

    draw_reconstruction_field!(
        fig,
        5,
        ps_ensemble_simulation,
        eof_result,
        member,
        :ps,
        modes,
        current_index;
        data_label="PS Values Pa",
        scope_index=scope_index,
        data_colormap=data_colormap,
        coastline_color=coastline_color,
        shading=shading,
        fontsize=fontsize,
        whole_figure_title="Anomalies of PS fields",
        display_mode=:anomaly,
        factor_for_data_transform=1
    )

    record(fig, filename, chosen_indices; framerate=framerate) do t
        current_index[] = t
    end


    return fig

end


function temporal_pattern_integration(eof_data::EOFEnsemble; mode=1, data_name="", size=Makie.autmatic, value_mode=:normal, transformation=nothing, legend_position=:rt)

    xmappings = Int[]

    yvals = Float64[]

    scopes = eof_data.scopes


    control_integration = Vector{Float64}(undef, length(scopes))

    xvals = eachindex(scopes)

    time_axis = eof_data.time



    function transform_data(data)
        if isnothing(transformation)
            return data
        else
            return standardize(transformation, data)
        end
    end

    for scope_index in eachindex(scopes)

        for member_id in keys(eof_data.ensemble)

            eof_temporalpattern = eof_data.ensemble[member_id][scope_index].temporal_modes[:, mode]


            integration_result = integrate(collect(1:length(eof_temporalpattern)), transform_data(eof_temporalpattern))


            push!(xmappings, scope_index)
            push!(yvals, integration_result)

        end

        control_integration[scope_index] = integrate(collect(1:length(scopes[scope_index])), transform_data(eof_data.piControl[scope_index].temporal_modes[:, mode]))
    end


    fig = Figure(size=size)

    ax_boxplot = Axis(fig[1, 1], title="Integration of $data_name", yminorticksvisible=true)

    boxplot!(ax_boxplot, xmappings, yvals, label="Ensemble Integration Boxplot")

    scatter!(ax_boxplot, xvals, control_integration; color=:red, label="piControl Simulation")

    axislegend(ax_boxplot, position=legend_position)
    seasonslice = 1:3:length(eof_data.scopes)

    ax_boxplot.xticks = (seasonslice, ["$(year(time_axis[scopes[i].start])) - $(year(time_axis[scopes[i].stop]))" for i in seasonslice])
    ax_boxplot.xticklabelrotation = π / 4
    ax_boxplot.xticklabelalign = (:right, :center)

    return fig


end

function draw_contour_mode(
    axis,
    eof_result,
    scope_index;
    single_color=nothing,
    hexbin_colormap=:algae,)

    contour_colors = distinguishable_colors(50)

    result = []

    if contour_mode == :spaghetti
        for (i, (_, eofres_array)) in enumerate(eof_result.ensemble)

            color = isnothing(single_color) ? contour_colors[i] : single_color
            contour = contour!(
                axis,
                lonmin .. lonmax,
                latmin .. latmax,
                @lift(eofres_array[$scope_index].spatial_modes[:, :, mode] * scale_value),
                levels=contour_levels,
                color=color,
                transformation=(; translation=(0, 0, 1100))
            )

            push!(result, contour)
        end
    elseif contour_mode == :hexbin

        picontrol_observable = @lift(eof_result.piControl[$current_scope_index].spatial_modes[:, :, mode] * scale_value)

        contour!(
            axis,
            lonmin .. lonmax,
            latmin .. latmax,
            picontrol_observable,
            levels=contour_levels,
            color=:red,
            transformation=(; translation=(0, 0, 1100))
        )


        contours_observable = @lift([contours(eof_result.lon, eof_result.lat, eofres_array[$current_scope_index].spatial_modes[:, :, mode], contour_levels) for (i, (_, eofres_array)) in enumerate(eof_result.ensemble)])

        vertices_observable = @lift(vcat(sample_along_line.(get_isocontour_vertices($contours_observable...); dx=0.01)...))

        hb = hexbin!(axis, vertices_observable, cellsize=2.0, colormap=hexbin_colormap, threshold=1)


    else
        ArgumentError("Use either :spaghetti or :hexbin")
    end

end


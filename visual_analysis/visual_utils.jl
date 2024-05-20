using GLMakie
using GeoMakie
using EmpiricalOrthogonalFunctions
using NCDatasets
using Dates
using BenchmarkTools
using Statistics
using StatsBase
using Colors

include("utils.jl")

# create geographical axes
function local_geoaxis_creation!(
    figure::Makie.Figure,
    lonlims::Tuple{<:Number,<:Number},
    latlims::Tuple{<:Number,<:Number};
    figure_row::Int=1,
    figure_col::Int=1,
    title="",
    dest_projection="+proj=merc",
)

    geoaxis = GeoAxis(
        figure[figure_row, figure_col];
        dest=dest_projection,
        limits=(lonlims, latlims),
        title=title,
    )

    return geoaxis
end


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
)::Makie.Figure

    nmodes = length(eof_data[get_member_id_string(1)][1].modes_variability)


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
        title = @lift("$(data.id) $(year(data.time[all_scopes[$current_scope_index].start]))-$(year(data.time[all_scopes[$current_scope_index].stop])) Mode $modenr Variability $(round(mean([resarray[$current_scope_index].modes_variability[modenr] for (_, resarray) in eof_data]), digits = 2)) %")
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

    contour_colors = distinguishable_colors(50)

    for (i, (_, eofres_array)) in enumerate(eof_data)
        for mode in 1:nmodes
            contour!(mode_axes[mode], lonmin .. lonmax, latmin .. latmax, @lift(eofres_array[$current_scope_index].spatial_modes[:, :, mode]), levels=[0], color=contour_colors[i], transformation=(; translation=(0, 0, 1100)))
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


function display_eof_modes!(
    ensemble_simulation::EnsembleSimulation,
    eof_data::Dict{String,Vector{EOFResult}},
    scopes,
    current_scope,
    axis_array,
    limits;
    contour_levels=[0],
    contours_for_slice=:,
    colormap=:viridis,
    coastline_color=:white,
    shading=Makie.automatic)

    lonmin, lonmax = extrema(ensemble_simulation.lons)
    latmin, latmax = extrema(ensemble_simulation.lats)


    nmodes = length(eof_data[get_member_id_string(1)][1].modes_variability)

    mean_axis = axis_array[end]

    dataset_mean_observable = @lift(dropdims(mean(get_mean_of_multiple_arrays([member.data[:, :, $current_scope] for member in data.members]...), dims=3), dims=3))

    surface!(
        mean_axis,
        lonmin .. lonmax,
        latmin .. latmax,
        dataset_mean_observable;
        shading=shading,
        colormap=colormap,
        colorrange=limits,
        overdraw=false,
        transformation=(; translation=(0, 0, -1000))
    )
    lines!(mean_axis, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

    mode_iterator = 1:nmodes

    contour_colors = distinguishable_colors(50)


    for mode in mode_iterator
        surface!(
            axis_array[mode],
            lonmin .. lonmax,
            latmin .. latmax,
            @lift(get_mean_of_multiple_arrays([res_array[$current_scope_index].spatial_modes[:, :, mode] for (_, res_array) in eof_data]...));
            shading=shading,
            colormap=colormap,
            colorrange=limits,
            overdraw=false,
            transformation=(; translation=(0, 0, -1 * limits[2]))
        )
        lines!(axis_array[mode], GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

        if mode in mode_iterator[contours_for_slice]
            for (i, (_, eofres_array)) in enumerate(eof_data)
                contour!(
                    axis_array[mode],
                    lonmin .. lonmax,
                    latmin .. latmax,
                    @lift(eofres_array[$current_scope_index].spatial_modes[:, :, mode]),
                    levels=contour_levels,
                    color=contour_colors[i],
                    transformation=(; translation=(0, 0, 1100))
                )
            end
        end

    end



end

function pick_correct_field(eof_ensemble_results::EOFEnsembleResult, field)
    if field == :ivt
        piControl_eofs = eof_ensemble_results.ivt_piControl
        ensemble_data = eof_ensemble_results.ivt_ensemble_eofs

    elseif field == :ps
        piControl_eofs = eof_ensemble_results.ps_piControl
        ensemble_data = eof_ensemble_results.ps_ensemble_eofs
    else
        ArgumentError("Use fields :ivt or :ps")
    end

    return ensemble_data, piControl_eofs
end

function compare_ensembles(
    nmodes,
    filename,
    ensembles_eof_tuples::Tuple{EnsembleSimulation,EOFEnsembleResult}...;
    field=:ivt,
    framerate::Int=30,
    colormap=:viridis,
    coastline_color=:white,
    shading=Makie.automatic,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing,
    fontsize::Int=12,
    whole_figure_title=nothing,
    contours_for_slice=:
)
    lonmin, lonmax = extrema(ensembles_eof_tuples[1][1].lons)
    latmin, latmax = extrema(ensembles_eof_tuples[1][1].lats)

    all_scopes = ensembles_eof_tuples[1][2].scopes 
    current_scope_index = Observable(1)

    fig = isnothing(resolution) ? Figure(fontsize=fontsize) : Figure(size=resolution, fontsize=fontsize)

    eof_extremas = Tuple{Float64,Float64}[]

    if !isnothing(whole_figure_title)
        Label(fig[0, 1:nmodes+2], whole_figure_title, fontsize=round(1.7 * fontsize))
    end

    for (_, eof_data) in ensembles_eof_tuples

        ensemble_eof_data, _ = pick_correct_field(eof_data, :ivt)
        
        append!(eof_extremas, [extrema(res.spatial_modes) for (_, eofresarray) in ensemble_eof_data for res in eofresarray])
    end

    limits = reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), eof_extremas)

    for (i, (ensemble, eof_data)) in enumerate(ensembles_eof_tuples)
        mode_axes = []

        current_row = i+1

        ensemble_eof_data, _ = pick_correct_field(eof_data, :ivt)

        for modenr in 1:nmodes
            title = @lift("$(ensemble.id) $(year(ensemble.time[all_scopes[$current_scope_index].start]))-$(year(ensemble.time[all_scopes[$current_scope_index].stop])) Mode $modenr Variability $(round(mean([resarray[$current_scope_index].modes_variability[modenr] for (_, resarray) in ensemble_eof_data]), digits = 2)) %")
            push!(mode_axes, local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=title, figure_row=current_row, figure_col=modenr))
        end
        push!(mode_axes, local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=@lift("$(ensemble.id) $(year(ensemble.time[all_scopes[$current_scope_index].start]))-$(year(ensemble.time[all_scopes[$current_scope_index].stop])) PIControl Simulation"), figure_row=current_row, figure_col=nmodes + 1))

        push!(mode_axes, local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=@lift("$(ensemble.id) $(year(ensemble.time[all_scopes[$current_scope_index].start]))-$(year(ensemble.time[all_scopes[$current_scope_index].stop])) Mean map"), figure_row=current_row, figure_col=nmodes + 2))

        display_eof_modes!(
            ensemble,
            eof_data,
            current_scope_index,
            mode_axes,
            limits;
            field=field,
            contour_levels=[0],
            contours_for_slice=contours_for_slice,
            colormap=colormap,
            coastline_color=coastline_color,
            shading=shading
        )

    end

    Colorbar(fig[1, 1:nmodes+2], limits=limits, colormap=colormap, vertical=false, label="EOF values")



    record(fig, filename, eachindex(all_scopes); framerate=framerate) do t
        current_scope_index[] = t
    end


    return fig



end

function display_eof_modes!(
    ensemble_simulation::EnsembleSimulation,
    eof_result::EOFEnsembleResult,
    current_scope_index,
    axis_array,
    limits;
    field=:ivt,
    contour_levels=[0],
    contours_for_slice=:,
    colormap=:viridis,
    coastline_color=:white,
    shading=Makie.automatic)

    lonmin, lonmax = extrema(ensemble_simulation.lons)
    latmin, latmax = extrema(ensemble_simulation.lats)

    ensemble_eofs, piControl_eofs = pick_correct_field(eof_result, field)


    nmodes = length(ensemble_eofs[get_member_id_string(1)][1].modes_variability)

    mean_axis = axis_array[end]
    picontrol_axis = axis_array[end-1]

    all_scopes = eof_result.scopes

    dataset_mean_observable = @lift(dropdims(mean(get_mean_of_multiple_arrays([member.data[:, :, all_scopes[$current_scope_index]] for member in ensemble_simulation.members]...), dims=3), dims=3))

    picontrol_observable = @lift(piControl_eofs[$current_scope_index].spatial_modes[:, :, 1])

    contour_colors = distinguishable_colors(50)

    surface!(
        picontrol_axis,
        lonmin .. lonmax,
        latmin .. latmax,
        picontrol_observable;
        shading=shading,
        colormap=colormap,
        colorrange=limits,
        overdraw=false,
        transformation=(; translation=(0, 0, -1000))
    )
    contour!(
        picontrol_axis,
        lonmin .. lonmax,
        latmin .. latmax,
        picontrol_observable,
        levels=contour_levels,
        color=:red,
        transformation=(; translation=(0, 0, 1100))
    )
    lines!(picontrol_axis, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

    surface!(
        mean_axis,
        lonmin .. lonmax,
        latmin .. latmax,
        dataset_mean_observable;
        shading=shading,
        colormap=colormap,
        colorrange=limits,
        overdraw=false,
        transformation=(; translation=(0, 0, -1000))
    )
    lines!(mean_axis, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

    mode_iterator = 1:nmodes


    for mode in mode_iterator
        surface!(
            axis_array[mode],
            lonmin .. lonmax,
            latmin .. latmax,
            @lift(get_mean_of_multiple_arrays([res_array[$current_scope_index].spatial_modes[:, :, mode] for (_, res_array) in ensemble_eofs]...));
            shading=shading,
            colormap=colormap,
            colorrange=limits,
            overdraw=false,
            transformation=(; translation=(0, 0, -1 * limits[2]))
        )
        lines!(axis_array[mode], GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

        if mode in mode_iterator[contours_for_slice]
            for (i, (_, eofres_array)) in enumerate(ensemble_eofs)
                contour!(
                    axis_array[mode],
                    lonmin .. lonmax,
                    latmin .. latmax,
                    @lift(eofres_array[$current_scope_index].spatial_modes[:, :, mode]),
                    levels=contour_levels,
                    color=contour_colors[i],
                    transformation=(; translation=(0, 0, 1100))
                )
            end
        end

    end



end



function generate_correlation_boxplot(ivt_eof_ensemble, ps_eof_ensemble, scopes, time_axis; compare_modes=1 => 1, data_name="", size=Makie.autmatic, control_data=nothing)

    available_members = intersect(keys(ps_eof_ensemble), keys(ivt_eof_ensemble))

    println("Available members in both datasets: $(length(available_members))")

    xmappings = Int[]

    yvals = Float64[]

    lagvals = Int[]

    lag = collect(-100:100)

    for scope_index in eachindex(scopes)

        for member_id in available_members

            ps_eof_temporalpattern = ps_eof_ensemble[member_id][scope_index].temporal_modes[:, compare_modes.second]
            ivt_eof_temporalpattern = ivt_eof_ensemble[member_id][scope_index].temporal_modes[:, compare_modes.first]

            # println("Length of signals: $(length(ps_eof_temporalpattern)) $(length(ivt_eof_temporalpattern))")
            # half_signal = trunc(Int, round(length(ivt_eof_temporalpattern) * 0.5))
            # lag = collect(-1*half_signal:half_signal)


            corsscor_result = crosscor(standardize(ZScoreTransform, ivt_eof_temporalpattern), standardize(ZScoreTransform, ps_eof_temporalpattern), lag)

            (maximal_correlation, index) = findmax(abs.(corsscor_result))

            push!(xmappings, scope_index)
            push!(yvals, maximal_correlation)
            push!(lagvals, lag[index])
        end
    end

    fig = Figure(size=size)

    ax_boxplot = Axis(fig[1, 1], title="Crosscorrelation of $data_name", yminorticksvisible=true)

    ax_lagval = Axis(fig[2, 1], title="Lags of Crosscorrelation of $data_name", yminorticksvisible=true)

    boxplot!(ax_boxplot, xmappings, yvals)
    boxplot!(ax_lagval, xmappings, lagvals)



    if !isnothing(control_data)

        xvals = eachindex(scopes)

        control_correlation = Vector{Float64}(undef, length(scopes))
        control_lag = Vector{Int}(undef, length(scopes))

        for i in eachindex(scopes)

            ivt_temp_signal = control_data[1][i]
            ps_temp_signal = control_data[2][i]

            corsscor_result = crosscor(standardize(ZScoreTransform, ivt_temp_signal), standardize(ZScoreTransform, ps_temp_signal), lag)

            (maximal_correlation, index) = findmax(abs.(corsscor_result))

            control_correlation[i] = maximal_correlation
            control_lag[i] = index
        end

        lines!(ax_boxplot, xvals, control_correlation; color=:red)
        lines!(ax_lagval, xvals, control_lag; color=:red)
    end




    seasonslice = 1:3:length(scopes)

    for axis in [ax_boxplot, ax_lagval]
        axis.xticks = (seasonslice, ["$(year(time_axis[scopes[i].start])) - $(year(time_axis[scopes[i].stop]))" for i in seasonslice])
        axis.xticklabelrotation = π / 4
        axis.xticklabelalign = (:right, :center)
        axis
    end



    return fig


end

function generate_correlation_boxplot(eof_result_data::EOFEnsembleResult, time_axis; compare_modes=1 => 1, data_name="", size=Makie.autmatic)

    available_members = intersect(keys(eof_result_data.ps_ensemble_eofs), keys(eof_result_data.ivt_ensemble_eofs))

    println("Available members in both datasets: $(length(available_members))")

    xmappings = Int[]

    yvals = Float64[]

    lagvals = Int[]

    control_correlation = Vector{Float64}(undef, length(eof_result_data.scopes))
    control_lag = Vector{Int}(undef, length(eof_result_data.scopes))

    xvals = eachindex(eof_result_data.scopes)

    for scope_index in eachindex(eof_result_data.scopes)

        lag_extent = Int(floor(length(eof_result_data.scopes[scope_index]) / 2))

        lag = collect(-1*lag_extent:lag_extent)

        for member_id in available_members

            ps_eof_temporalpattern = eof_result_data.ps_ensemble_eofs[member_id][scope_index].temporal_modes[:, compare_modes.second]
            ivt_eof_temporalpattern = eof_result_data.ivt_ensemble_eofs[member_id][scope_index].temporal_modes[:, compare_modes.first]

            # println("Length of signals: $(length(ps_eof_temporalpattern)) $(length(ivt_eof_temporalpattern))")
            # half_signal = trunc(Int, round(length(ivt_eof_temporalpattern) * 0.5))
            # lag = collect(-1*half_signal:half_signal)


            corsscor_result = crosscor(standardize(ZScoreTransform, ivt_eof_temporalpattern), standardize(ZScoreTransform, ps_eof_temporalpattern), lag)

            (maximal_correlation, index) = findmax(abs.(corsscor_result))

            push!(xmappings, scope_index)
            push!(yvals, maximal_correlation)
            push!(lagvals, lag[index])
        end



        corsscor_result = crosscor(standardize(ZScoreTransform, eof_result_data.ivt_piControl[scope_index].temporal_modes[:, compare_modes.first]), standardize(ZScoreTransform, eof_result_data.ps_piControl[scope_index].temporal_modes[:, compare_modes.second]), lag)

        (maximal_correlation, index) = findmax(abs.(corsscor_result))

        control_correlation[scope_index] = maximal_correlation
        control_lag[scope_index] = lag[index]
    end

    fig = Figure(size=size)

    ax_boxplot = Axis(fig[1, 1], title="Crosscorrelation of $data_name", yminorticksvisible=true)

    ax_lagval = Axis(fig[2, 1], title="Lags of Crosscorrelation of $data_name", yminorticksvisible=true)
    scatter!(ax_lagval, xvals, control_lag; color=:red, label="piControl simulation")

    boxplot!(ax_boxplot, xmappings, yvals)
    boxplot!(ax_lagval, xmappings, lagvals)

    scatter!(ax_boxplot, xvals, control_correlation; color=:red, label="piControl simulation")

    axislegend(ax_boxplot)
    axislegend(ax_lagval)
    seasonslice = 1:3:length(eof_result_data.scopes)

    for axis in [ax_boxplot, ax_lagval]
        axis.xticks = (seasonslice, ["$(year(time_axis[eof_result_data.scopes[i].start])) - $(year(time_axis[eof_result_data.scopes[i].stop]))" for i in seasonslice])
        axis.xticklabelrotation = π / 4
        axis.xticklabelalign = (:right, :center)
        axis
    end



    return fig


end


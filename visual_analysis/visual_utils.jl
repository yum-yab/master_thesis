using GLMakie
using GeoMakie
using EmpiricalOrthogonalFunctions
using NCDatasets
using Dates
using BenchmarkTools
using Statistics
using StatsBase
using Colors

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
        pc_axis[dataset.name] = Axis(fig[2*i, 1:nmodes+1], limits = ((nothing, nothing), (-1, 1)))
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


    pc_axis = Axis(fig[3, 1:nmodes+1], limits = ((nothing, nothing), (-1, 1)))

    lonmin, lonmax = extrema(data.lons)
    latmin, latmax = extrema(data.lats)

    eof_extremas = Tuple{Float64,Float64}[]

    pcs_extremas = Tuple{Float64,Float64}[]

    for (_, eofres_array) in eof_data
        append!(eof_extremas, [extrema(res.spatial_modes) for res in eofres_array])
    end

    spatial_min_val, spatial_max_val = reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), eof_extremas)
    # spatial_min_val = -1 * spatial_max_val
    println("Data extremas: $spatial_min_val $spatial_max_val")


    mode_axes = []

    for modenr in 1:nmodes
        title = @lift("$(data.id) $(year(data.time[all_scopes[$current_scope_index].start]))-$(year(data.time[all_scopes[$current_scope_index].stop])) Mode $modenr Variability $(round(mean([resarray[$current_scope_index].modes_variability[modenr] for (_, resarray) in eof_data]), digits = 2)) %")
        push!(mode_axes, local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=title, figure_row=1, figure_col=modenr))
    end

    mean_axis = local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=@lift("$(data.id) $(year(data.time[all_scopes[$current_scope_index].start]))-$(year(data.time[all_scopes[$current_scope_index].stop])) Mean map"), figure_row=1, figure_col=nmodes+1)

    dataset_mean_observable = @lift(dropdims(mean(get_mean_of_multiple_arrays([member.data[:, :, all_scopes[$current_scope_index]] for member in data.members]...), dims = 3), dims = 3))

    surface!(mean_axis, lonmin .. lonmax, latmin .. latmax, dataset_mean_observable; shading=shading, colormap=colormap, colorrange=(spatial_min_val, spatial_max_val), alpha=surface_alpha_value, overdraw=false, transformation=(; translation=(0, 0, -1000)))
    lines!(mean_axis, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))


    for mode in 1:nmodes
        surface!(mode_axes[mode], lonmin .. lonmax, latmin .. latmax, @lift(get_mean_of_multiple_arrays([res_array[$current_scope_index].spatial_modes[:, :, mode] for (_, res_array) in eof_data]...)); shading=shading, colormap=colormap, colorrange=(spatial_min_val, spatial_max_val), alpha=surface_alpha_value, overdraw=false, transformation=(; translation=(0, 0, -1 * spatial_max_val)))
        lines!(mode_axes[mode], GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))
        positions = @lift(collect(zip(eachindex(all_scopes[$current_scope_index]), standardize(UnitRangeTransform, get_mean_of_multiple_arrays([res_array[$current_scope_index].temporal_modes[:, mode] for (_, res_array) in eof_data]...); unit=false))))
        lines!(pc_axis, positions; label="Mode $mode")
    end
    axislegend(pc_axis)

    contour_colors = distinguishable_colors(50)

    for (i, (_, eofres_array)) in enumerate(eof_data)
        for mode in 1:nmodes
            contour!(mode_axes[mode], lonmin .. lonmax, latmin .. latmax, @lift(eofres_array[$current_scope_index].spatial_modes[:, :, mode]), levels=[0], color= contour_colors[i], transformation=(; translation=(0, 0, 1100)))
        end
    end



    Colorbar(fig[2, 1:nmodes+1], limits=(spatial_min_val, spatial_max_val), colormap=colormap, vertical=false, label="EOF values")

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

function generate_correlation_boxplot(ivt_eof_ensemble, ps_eof_ensemble, scopes, time_axis; compare_modes = 1 => 1, title = "", size = Makie.autmatic)

    available_members = intersect(keys(ps_eof_ensemble), keys(ivt_eof_ensemble))

    println("Available members in both datasets: $(length(available_members))")

    xmappings = Int[]

    yvals = Float64[]

    lagvals = Int[]

    for scope in eachindex(scopes)

        for member_id in available_members

            ps_eof_temporalpattern = ps_eof_ensemble[member_id][scope].temporal_modes[:, compare_modes.second]
            ivt_eof_temporalpattern = ivt_eof_ensemble[member_id][scope].temporal_modes[:, compare_modes.first]
            
            println("Length of signals: $(length(ps_eof_temporalpattern)) $(length(ivt_eof_temporalpattern))")
            half_signal = trunc(Int, round(length(ivt_eof_temporalpattern) * 0.5))
            lag = collect(-1 *half_signal:half_signal)
            lag = collect(-100:100)

            corsscor_result = crosscor(standardize(ZScoreTransform, ivt_eof_temporalpattern), standardize(ZScoreTransform, ps_eof_temporalpattern), lag)

            (maximal_correlation, index) = findmax(abs.(corsscor_result))

            push!(xmappings, scope)
            push!(yvals, maximal_correlation)
            push!(lagvals, lag[index])
        end
    end

    println("Index of crosscor: $index_set")

    fig = Figure(size = size)

    ax_boxplot = Axis(fig[1,1], title = title)

    ax_lagval = Axis(fig[2,1], title = title)

    boxplot!(ax_boxplot, xmappings, yvals)
    boxplot!(ax_lagval, xmappings, lagvals)
    
    seasonslice = 1:3:length(scopes)

    for axis in [ax_boxplot, ax_lagval]
        axis.xticks = (seasonslice, ["$(year(time_axis[scopes[i].start])) - $(year(time_axis[scopes[i].stop]))" for i in seasonslice])
        axis.xticklabelrotation = π / 4
        axis.xticklabelalign = (:right, :center)
    end

    

    return fig
    
    
end
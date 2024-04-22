using GLMakie
using GeoMakie
using EmpiricalOrthogonalFunctions
using NCDatasets 
using Dates
using BenchmarkTools
using Statistics

# create geographical axes
function local_geoaxis_creation!(
    figure::Makie.Figure,
    lonlims::Tuple{<:Number, <:Number}, 
    latlims::Tuple{<:Number, <:Number}; 
    figure_row::Int = 1,
    figure_col::Int = 1,
    title = "",
    dest_projection = "+proj=merc",
    )

    geoaxis = GeoMakie.GeoAxis(
        figure[figure_row, figure_col];
        dest = dest_projection,
        limits=(lonlims, latlims),
        title = title
    )

    return geoaxis
end


function build_figure(
    data,
    lon_bounds::Tuple{<:Number, <:Number},
    lat_bounds::Tuple{<:Number, <:Number};
    title::String = "",
    colormap = :viridis,
    shading = Makie.automatic,
    resolution::Union{Nothing, Tuple{Int, Int}} = nothing
    )::Makie.Figure

    fig = isnothing(resolution) ?  Figure(fontsize=12) : Figure(size = resolution, fontsize=12)

    minval, maxval = extrema(data)

    ga = local_geoaxis_creation!(fig, lon_bounds, lat_bounds; title=title)
    surface!(ga, lon_bounds[1]..lon_bounds[2], lat_bounds[1]..lat_bounds[2], data; shading = shading, colormap = colormap)
    lines!(ga, GeoMakie.coastlines(); color = :white, transformation = (; translation = (0, 0, 1000)))
    Colorbar(fig[1, 2], limits = (minval, maxval), colormap = colormap)
    return fig
end

function animate_timeline(
    data::TimelineData,
    filename::String;
    framerate::Int = 30,
    colormap = :viridis,
    coastline_color = :white,
    shading = Makie.automatic,
    resolution::Union{Nothing, Tuple{Int, Int}} = nothing
    )::Makie.Figure

    fig = isnothing(resolution) ?  Figure(fontsize=12) : Figure(size = resolution, fontsize=12)

    time_index = Observable(1)

    axis = Dict()

    lon_bounds = extrema(data.lons)
    lat_bounds = extrema(data.lats)

    all_extrema = extrema.([d.data for d in data.datasets])

    min_val, max_val = reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), all_extrema)
    Colorbar(fig[2, 1:2], limits = (min_val, max_val), colormap = colormap, vertical = false,  label = "IVT (norm) in kg s^-1 m^-1")

    for (i, dataset) in enumerate(data.datasets)
        axis[dataset.name] = local_geoaxis_creation!(fig, lon_bounds, lat_bounds; title=@lift("$(dataset.name) $(Dates.format(data.time[$time_index], "mm/yyyy"))"), figure_col=i)
    end
    
    for dataset in data.datasets

        array_slice = @lift(dataset.data[:, :, $time_index])
        
        surface!(axis[dataset.name], lon_bounds[1]..lon_bounds[2], lat_bounds[1]..lat_bounds[2], array_slice; shading = shading, colormap = colormap, colorrange = (min_val, max_val))
        lines!(axis[dataset.name], GeoMakie.coastlines(); color = coastline_color, transformation = (; translation = (0, 0, 1000)))
    end

    timestamps = eachindex(data.time)
    record(fig, filename, timestamps; framerate = framerate) do t
        time_index[] = t
    end

    
    return fig
end

function show_eof_modes_of_timeline(
    data::TimelineData,
    nmodes::Int,
    n_seasons::Int,
    filename::String;
    eof_center::Bool = true,
    framerate::Int = 30,
    colormap = :viridis,
    coastline_color = :white,
    colorrange::Union{Nothing, Tuple{<:Real, <:Real}} = nothing, 
    shading = Makie.automatic,
    resolution::Union{Nothing, Tuple{Int, Int}} = nothing
    )::Makie.Figure

    fig = isnothing(resolution) ?  Figure(fontsize=12) : Figure(size = resolution, fontsize=12)

    all_scopes = get_sliding_time_scopes_by_threshold(data.time, n_seasons)

    current_scope = Observable(all_scopes[1])

    axis = Dict()

    lon_bounds = extrema(data.lons)
    lat_bounds = extrema(data.lats)

    Colorbar(fig[3, 1:nmodes], limits = colorrange, colormap = colormap, vertical = false,  label = "EOF values")

    for (i, dataset) in enumerate(data.datasets)
        all_modes_axis = []

        for j in 1:nmodes
            title = @lift("$(dataset.name) $(year(data.time[$current_scope.start]))-$(year(data.time[$current_scope.stop])) Mode $j")
            push!(all_modes_axis, local_geoaxis_creation!(fig, lon_bounds, lat_bounds; title=title, figure_row=i, figure_col=j))
        end
        
        axis[dataset.name] = all_modes_axis
    end
    
    for dataset in data.datasets

        eof_observable = @lift(get_eof_of_datachunk(dataset.data[:, :, $current_scope]; nmodes = nmodes))
        
        for mode in 1:nmodes
            surface!(axis[dataset.name][mode], lon_bounds[1]..lon_bounds[2], lat_bounds[1]..lat_bounds[2], @lift($eof_observable[1][:, :, mode]); shading = shading, colormap = colormap)
            lines!(axis[dataset.name][mode], GeoMakie.coastlines(); color = coastline_color, transformation = (; translation = (0, 0, 1000)))
        end
    end

    record(fig, filename, all_scopes; framerate = framerate) do t
        current_scope[] = t
    end

    
    return fig
end
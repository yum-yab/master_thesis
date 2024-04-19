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

function show_eof_modes_of_timeline(
    data::TimelineData,
    nmodes::Int,
    n_seasons::Int,
    filename::String;
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

    for (i, dataset) in enumerate(data.datasets)
        all_modes_axis = []

        for j in 1:nmodes
            title = @lift("$(dataset.name) $(year(data.time[$current_scope.start]))-$(year(data.time[$current_scope.stop])) Mode $j")
            push!(all_modes_axis, local_geoaxis_creation!(fig, lon_bounds, lat_bounds; title=title, figure_row=i, figure_col=j))
        end
        
        axis[dataset.name] = all_modes_axis
    end
    
    for dataset in data.datasets

        datachunk = @lift(dataset.data[:, :, $current_scope])

        spatial_modes, temporal_modes, variability_percentage = get_eof_of_datachunk(to_value(datachunk); nmodes = nmodes)

        for mode in 1:nmodes
            surface!(axis[dataset.name][mode], lon_bounds[1]..lon_bounds[2], lat_bounds[1]..lat_bounds[2], spatial_modes[:, :, mode]; shading = shading, colormap = colormap)
            lines!(axis[dataset.name][mode], GeoMakie.coastlines(); color = coastline_color, transformation = (; translation = (0, 0, 1000)))
        end
    end

    record(fig, filename, all_scopes; framerate = framerate) do t
        current_scope[] = t
    end

    
    return fig
end
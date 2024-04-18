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
    source_projection = Makie.automatic,
    dest_projection = "+proj=merc",
    )

    geoaxis = GeoMakie.GeoAxis(
        figure[figure_row, figure_col];
        dest = dest_projection,
        source = source_projection,
        limits=(lonlims, latlims),
        title = title
    )

    return geoaxis
end
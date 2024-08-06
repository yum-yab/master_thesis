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

# create geographical axes
function local_geoaxis_creation!(
    figure::Makie.Figure,
    lonlims::Tuple{<:Number,<:Number},
    latlims::Tuple{<:Number,<:Number};
    figure_row=1,
    figure_col=1,
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






function compare_ensembles(
    nmodes,
    filename,
    ensembles_eof_tuples::Tuple{EOFEnsemble,Float64,Any}...;
    framerate::Int=30,
    coastline_color=:white,
    shading=Makie.automatic,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing,
    fontsize::Int=12,
    whole_figure_title=nothing,
    contours_for_slice=:,
    contour_mode=:spaghetti,
    hexbin_colormap=:amp,
    testmode=false,
    horizontal=true
)

    all_scopes = ensembles_eof_tuples[1][1].scopes
    time_axis = ensembles_eof_tuples[1][1].time
    current_scope_index = Observable(1)

    current_scope = @lift(all_scopes[$current_scope_index])

    fig = isnothing(resolution) ? Figure(fontsize=fontsize) : Figure(size=resolution, fontsize=fontsize)

    if !isnothing(whole_figure_title)
        date_str_obs = @lift("$(year(time_axis[$current_scope.start]))-$(year(time_axis[$current_scope.stop]))")
        Label(fig[0, 1:nmodes+2], @lift("$whole_figure_title $($date_str_obs)"), fontsize=round(1.7 * fontsize))
    end


    for (i, (eof_data, data_scaling, colormap)) in enumerate(ensembles_eof_tuples)

        current_row = i

        display_eof_modes!(
            fig,
            eof_data,
            current_scope_index,
            current_row,
            nmodes;
            contour_levels=[0],
            contours_for_slice=contours_for_slice,
            colormap=colormap,
            coastline_color=coastline_color,
            shading=shading,
            contour_mode=contour_mode,
            hexbin_colormap=hexbin_colormap,
            scale_value=data_scaling,
            use_rows=horizontal
        )

    end


    if !testmode
        record(fig, filename, eachindex(all_scopes); framerate=framerate) do t
            current_scope_index[] = t
        end
    end

    return fig



end


function display_eof_modes!(
    layout,
    eof_result::EOFEnsemble,
    current_scope_index,
    nmodes;
    contour_levels=[0],
    contour_mode=:spaghetti,
    hexbin_colormap=:amp,
    colormap=:viridis,
    coastline_color=:white,
    shading=Makie.automatic,
    vertical=true,
    sample_distance = 0.1,
    axis_size = 500
)
    lonmin, lonmax = extrema(eof_result.lon)
    latmin, latmax = extrema(eof_result.lat)

 
    picontrol_observable = @lift(eof_result.piControl[$current_scope_index].spatial_modes[:, :, 1])

    contour_colors = distinguishable_colors(50)

    eof_extremas = [extrema(res.spatial_modes) for (_, eofresarray) in eof_result.ensemble for res in eofresarray]

    _, maxval = reduce((a, b) -> (min(a[1], b[1]), max(a[2], b[2])), eof_extremas)

    limits = (-1 * maxval, maxval)

    mode_iterator = 1:nmodes

    function get_position(mode_id)
        if vertical
            return (mode_id, 1)
        else
            return (1, mode_id)
        end
    end

    for mode in mode_iterator

        row, column = get_position(mode)

        mode_layout = layout[row, column] = GridLayout()

        # title = @lift("$(ensemble_simulation.id) $(year(ensemble_simulation.time[all_scopes[$current_scope_index].start]))-$(year(ensemble_simulation.time[all_scopes[$current_scope_index].stop])) Mode $mode Variability $(round(mean([get_modes_variability(resarray[$current_scope_index])[mode] for (_, resarray) in eof_result.ensemble]), digits = 2) * 100) %")
        title = "Mode $mode"
        axis = GeoAxis(
            mode_layout[1, 1:2];
            dest="+proj=merc",
            limits=((lonmin, lonmax), (latmin, latmax)),
            title=title,
            width=axis_size, height=axis_size,
            # tellwidth=false, tellheight=false
        )
        # axis = local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=title, figure_row=row, figure_col=column)

        surf = surface!(
            axis,
            lonmin .. lonmax,
            latmin .. latmax,
            @lift(get_mean_of_multiple_arrays([res_array[$current_scope_index].spatial_modes[:, :, mode] for (_, res_array) in eof_result.ensemble]...));
            shading=shading,
            colormap=colormap,
            colorrange=limits,
            overdraw=false,
            transformation=(; translation=(0, 0, -1 * limits[2]))
        )
        lines!(axis, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))


        surface_toggle = Toggle(mode_layout[2, 1], active=true)
        ensemble_vis_toggle = Toggle(mode_layout[2, 2], active=true, buttoncolor = :red)

        connect!(surf.visible, surface_toggle.active)

        if contour_mode == :spaghetti

            seperate_matrices_observables = [@lift(eofres_array[$current_scope_index].spatial_modes[:, :, $mode]) for (_, eofres_array) in eof_result.ensemble]

            spaghetti_contours = [contour!(
                axis,
                lonmin .. lonmax,
                latmin .. latmax,
                matrix,
                levels=contour_levels,
                color=contour_colors[i],
                transformation=(; translation=(0, 0, 1100))
            ) for (i, matrix) in enumerate(seperate_matrices_observables)]

            for c in spaghetti_contours
                connect!(c.visible, ensemble_vis_toggle.active)
            end
            
        elseif contour_mode == :hexbin

            picontrol_observable = @lift(eof_result.piControl[$current_scope_index].spatial_modes[:, :, mode])

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

            ensemble_vertices_obs = @lift([vcat(sample_along_line.(lines; dx=sample_distance)...) for lines in get_isocontour_lines($contours_observable...)])

            ehb = ensemblehexbin!(axis, ensemble_vertices_obs, cellsize=2.0, colormap=hexbin_colormap, threshold=1)

            connect!(ehb.visible, ensemble_vis_toggle.active)


        else
            ArgumentError("Use either :spaghetti or :hexbin")
        end

        if !vertical
            colsize!(layout, mode, Fixed(axis_size))
        end

    end

    # for i in 1:nmodes
    #     colsize!(layout, mode, Fixed(axis_size))
    # end
    

    hexbin_cb_row, hexbin_cb_col = get_position(nmodes + 1)

    Colorbar(fig[hexbin_cb_row, hexbin_cb_col], colormap=hexbin_colormap, label="No. Members Contour Line hitting bin", colorrange=(0, length(eof_result.ensemble)), vertical=!vertical)



    cb_row, cb_col = get_position(nmodes + 2)

    Colorbar(fig[cb_row, cb_col], limits=limits, colormap=colormap, vertical=!vertical, label=get_var_unit_string(eof_result))


end

function generate_correlation_boxplot(
    first_eof_data::EOFEnsemble,
    second_eof_data::EOFEnsemble;
    compare_modes=1 => 1,
    data_name="",
    size=Makie.autmatic,
    value_mode=:normal,
    transformation=nothing,
    show_first_pattern=false,
    colormaps=(:vik100, :vik100),
    use_crosscor=true,
    fontsize=18)

    available_members = intersect(keys(first_eof_data.ensemble), keys(second_eof_data.ensemble))

    xmappings = Int[]

    yvals = Float64[]

    lagvals = Int[]

    scopes = first_eof_data.scopes


    control_correlation = Vector{Float64}(undef, length(scopes))
    control_lag = Vector{Int}(undef, length(scopes))

    xvals = eachindex(scopes)

    time_axis = first_eof_data.time

    member_colors = repeat(distinguishable_colors(50), length(scopes))

    function transform_data(data)
        if isnothing(transformation)
            return data
        else
            return standardize(transformation, data)
        end
    end

    lag_extent = 4

    lag = use_crosscor ? collect(-1*lag_extent:lag_extent) : [0]

    for scope_index in eachindex(scopes)

        # lag_extent = Int(floor(length(scopes[scope_index]) / 2))


        for member_id in available_members

            first_eof_temporalpattern = first_eof_data.ensemble[member_id][scope_index].temporal_modes[:, compare_modes.first]
            second_eof_temporalpattern = second_eof_data.ensemble[member_id][scope_index].temporal_modes[:, compare_modes.second]


            if length(first_eof_temporalpattern) != (length(second_eof_temporalpattern))
                println("Scope $(scope_index): Size fisrt: $(length(first_eof_temporalpattern))\tSize Second: $(length(second_eof_temporalpattern))")
                continue
            end

            # println("Length of signals: $(length(ps_eof_temporalpattern)) $(length(ivt_eof_temporalpattern))")
            # half_signal = trunc(Int, round(length(ivt_eof_temporalpattern) * 0.5))
            # lag = collect(-1*half_signal:half_signal)

            crosscor_result = crosscor(transform_data(first_eof_temporalpattern), transform_data(second_eof_temporalpattern), lag)

            (maximal_correlation, index) = findmax(abs.(crosscor_result))

            push!(xmappings, scope_index)

            if value_mode == :normal
                push!(yvals, crosscor_result[index])
            elseif value_mode == :absolute
                push!(yvals, maximal_correlation)
            else
                ArgumentError("Please use :normal or :absolute as value_mode")
            end
            push!(lagvals, lag[index])
        end



        crosscor_piControl_result = crosscor(transform_data(first_eof_data.piControl[scope_index].temporal_modes[:, compare_modes.first]), transform_data(second_eof_data.piControl[scope_index].temporal_modes[:, compare_modes.second]), lag)

        (maximal_correlation, index) = findmax(abs.(crosscor_piControl_result))

        if value_mode == :normal
            control_correlation[scope_index] = crosscor_piControl_result[index]
        elseif value_mode == :absolute
            control_correlation[scope_index] = maximal_correlation
        else
            ArgumentError("Please use :normal or :absolute as value_mode")
        end
        control_lag[scope_index] = lag[index]
    end

    if mean(yvals) > 0
        legend_position = :rb
    else
        legend_position = :rt
    end

    fig = Figure(size=size, fontsize=fontsize)

    boxplot_cols = use_crosscor ? UnitRange(1, 2) : UnitRange(1:4)

    value_type = value_mode == :absolute ? "Absolute" : ""

    cor_type = use_crosscor ? "Maximum $value_type Crosscorrelation" : "$value_type Correlation"

    ax_boxplot = Axis(fig[1, boxplot_cols], title="$cor_type of $data_name", yminorticksvisible=true, limits=(nothing, (-1, 1)), ylabel="Pearson correlation coefficient")

    if use_crosscor
        ax_lagval = Axis(fig[2, boxplot_cols], title="Lags of $cor_type of $data_name", yminorticksvisible=true, limits=(nothing, extrema(lag)))
        boxplot!(ax_lagval, xmappings, lagvals, outliercolor=member_colors, label="Lag Boxplot")
        scatter!(ax_lagval, xvals, control_lag; color=:red, label="piControl simulation")
        axislegend(ax_lagval, position=legend_position)
    end

    boxplot!(ax_boxplot, xmappings, yvals, outliercolor=member_colors, label="Ensemble Correlation Boxplot")
    scatter!(ax_boxplot, xvals, control_correlation; color=:red, label="piControl simulation")

    axislegend(ax_boxplot, position=legend_position)
    seasonslice = 1:5:length(first_eof_data.scopes)

    boxplot_axes = use_crosscor ? [ax_boxplot, ax_lagval] : [ax_boxplot]

    for axis in boxplot_axes
        axis.xticks = (seasonslice, ["$(year(time_axis[scopes[i].start])) - $(year(time_axis[scopes[i].stop]))" for i in seasonslice])
        axis.xticklabelrotation = π / 4
        axis.xticklabelalign = (:right, :center)
        axis
    end

    if show_first_pattern

        for (i, eof_data) in enumerate([first_eof_data, second_eof_data])

            index = 1

            lonmin, lonmax = extrema(eof_data.lon)
            latmin, latmax = extrema(eof_data.lat)

            time_axis = eof_data.time

            current_scope = eof_data.scopes[index]

            x, y = use_crosscor ? (i, 3) : (2, (2i - 1))


            colormap = colormaps[i]

            mode = compare_modes[i]

            title = "Mean Pattern of Mode $(mode) of $(uppercase(eof_data.variable_id)) Var. $(round(mean([get_modes_variability(resarray[1])[mode] for (_, resarray) in eof_data.ensemble]) * 100, digits = 2)) % $(year(time_axis[current_scope.start]))-$(year(time_axis[current_scope.stop]))"

            geoax = local_geoaxis_creation!(fig, (lonmin, lonmax), (latmin, latmax); title=title, figure_row=x, figure_col=y)


            mean_pattern = get_mean_of_multiple_arrays([res_array[index].spatial_modes[:, :, mode] for (_, res_array) in eof_data.ensemble]...)

            max_mode_1 = maximum(get_mean_of_multiple_arrays([res_array[index].spatial_modes[:, :, 1] for (_, res_array) in eof_data.ensemble]...))

            limits = (-max_mode_1, max_mode_1)
            surface!(
                geoax,
                lonmin .. lonmax,
                latmin .. latmax,
                mean_pattern;
                shading=NoShading,
                colormap=colormap,
                colorrange=limits,
                overdraw=false,
                transformation=(; translation=(0, 0, -1 * limits[2]))
            )
            lines!(geoax, GeoMakie.coastlines(); color=:black, transformation=(; translation=(0, 0, 1000)))

            cb_x, cb_y = use_crosscor ? (i, 4) : (2, 1 + (2i - 1))


            Colorbar(fig[cb_x, cb_y], limits=limits, colormap=colormap, vertical=true, label=get_var_unit_string(eof_data))
        end

    end


    return fig


end


function interactive_mode_analysis(
    nmodes,
    ensembles::EOFEnsemble...;
    coastline_color=:black,
    shading=Makie.automatic,
    resolution::Union{Nothing,Tuple{Int,Int}}=nothing,
    fontsize::Int=12,
    whole_figure_title=nothing,
    contour_mode=:spaghetti,
    hexbin_colormap=:algae,
    vertical=true,
    sample_distance = 0.1,
    axis_size=500
)

    fig = Figure()

    slider_grid = SliderGrid(
        fig[2, 1],
        (label="Scope", range=1:1:length(ensembles[1].scopes), format="{:%03d}", startvalue=1),
        width=1000,
        # tellheight=false,
        # tellwidth=false
    )

    current_scope_index = lift(slider_grid.sliders[1].value) do v
        return v
    end

    modes_layout = fig[1, 1] = GridLayout()

    for (i, eof_res) in enumerate(ensembles)

        r, c = vertical ? (1, i) : (i, 1)

        _, _, unit, cmap = get_correct_var_display(eof_res.variable_id)

        display_eof_modes!(
            modes_layout[r, c],
            eof_res,
            current_scope_index,
            nmodes;
            contour_levels=[0],
            contour_mode=contour_mode,
            hexbin_colormap=hexbin_colormap,
            colormap=cmap,
            coastline_color=coastline_color,
            shading=shading,
            vertical=vertical,
            sample_distance=sample_distance,
            axis_size=axis_size
        )

        if vertical
            colsize!(modes_layout, i, Fixed(axis_size))
        else
            rowsize!(modes_layout, i, Fixed(axis_size))
        end
        # colsize!(modes_layout, i, Fixed(axis_size))

    end

    


    resize_to_layout!(fig)
    return fig
end

function compare_ensemble_modes_variability(ensembles::Tuple{EOFEnsemble,String}...; resolution=(1080, 1920), fontsize=12, al_position=:rb)



    fig = Figure(size=resolution, fontsize=fontsize)

    all_axis = []


    for (i, (eof_ensemble, data_name)) in enumerate(ensembles)

        scopes = eof_ensemble.scopes

        time_axis = eof_ensemble.time

        nmodes = length(get_modes_variability(eof_ensemble.ensemble[get_member_id_string(11)][1]))

        xmappings = Int[]

        ymappings = [Float64[] for _ in 1:nmodes]

        piControly = [Float64[] for _ in 1:nmodes]

        for scope_index in eachindex(scopes)

            for (_, eof_results) in eof_ensemble.ensemble

                variability = get_modes_variability(eof_results[scope_index]) * 100

                push!(xmappings, scope_index)

                for mode in 1:nmodes
                    push!(ymappings[mode], variability[mode])
                end

            end

            piControl_var = get_modes_variability(eof_ensemble.piControl[scope_index]) * 100

            for mode in 1:nmodes
                push!(piControly[mode], piControl_var[mode])
            end
        end

        axis = Axis(fig[i, 1], title="In Modes Encoded Variability of $data_name", ytickformat="{:.1f} %")

        push!(all_axis, axis)

        # for mode in 1:nmodes
        #     println(size(piControly[mode]))
        # end

        for mode in 1:nmodes
            boxplot!(axis, xmappings, ymappings[mode], label="Mode $mode")
        end
        for mode in 1:nmodes
            lines!(axis, eachindex(scopes), piControly[mode], color=:red, label="piControl")
        end



        seasonslice = 1:3:length(scopes)

        axis.xticks = (seasonslice, ["$(year(time_axis[scopes[i].start])) - $(year(time_axis[scopes[i].stop]))" for i in seasonslice])
        axis.xticklabelrotation = π / 4
        axis.xticklabelalign = (:right, :center)


        axislegend(axis, position=al_position, unique=true)
    end


    # Label(fig[0, 1], "In Modes Encoded Variability of $data_name")

    linkyaxes!(all_axis...)



    return fig
end


function eof_data_correlation_maps!(
    layout,
    eof_result::EOFEnsemble,
    piControl::EnsembleSimulation,
    all_correlation_data::Vector{Array{Float64,4}};
    hexbin_colormap=:amp,
    colormap=:viridis,
    coastline_color=:white,
    shading=Makie.automatic,
    sample_distance=0.1
)

    nmodes = length(eof_result.ensemble[get_member_id_string(1)][1].singularvals)
    nmember = length(keys(eof_result.ensemble))

    lonmin, lonmax = extrema(eof_result.lon)
    latmin, latmax = extrema(eof_result.lat)


    all_scopes = eof_result.scopes

    contour_colors = distinguishable_colors(50)


    limits = (-1, 1)

    data_label = "Pearson Correlation"

    slider_grid = SliderGrid(
        layout[2, 1],
        (label="Contour Level", range=-1:0.1:1, format="{:.1f}", startvalue=0),
        (label="Scope", range=1:1:length(all_scopes), format="{:%03d}", startvalue=1),
        (label="Mode", range=1:1:nmodes, format="{:%1d}", startvalue=1),
        width=1000,
        # tellheight=false,
        # tellwidth=false
    )

    contour_levels = lift(slider_grid.sliders[1].value) do v
        return [v]
    end

    current_scope_index = lift(slider_grid.sliders[2].value) do v
        return v
    end

    mode = lift(slider_grid.sliders[3].value) do v
        return v
    end

    # all_correlation_data = @lift(generate_all_correlation_results(eof_result, original_data, $mode))





    mode_layout = layout[1, 1] = GridLayout()







    # title = @lift("$(ensemble_simulation.id) $(year(ensemble_simulation.time[all_scopes[$current_scope_index].start]))-$(year(ensemble_simulation.time[all_scopes[$current_scope_index].stop])) Mode $mode Variability $(round(mean([get_modes_variability(resarray[$current_scope_index])[mode] for (_, resarray) in eof_result.ensemble]), digits = 2) * 100) %")
    # title = @lift("Mode $mode Mean Variability: $(round(mean([get_modes_variability(resarray[$current_scope_index])[$mode] for (_, resarray) in eof_result.ensemble]) * 100, digits = 2)) %")
    title = @lift("Mode $($mode) $(year(eof_result.time[all_scopes[$current_scope_index].start]))-$(year(eof_result.time[all_scopes[$current_scope_index].stop]))")

    axis = GeoAxis(mode_layout[1, 1:3]; dest="+proj=merc", limits=((lonmin, lonmax), (latmin, latmax)), title=title,
        width=1000, height=1000,
        # tellwidth=false, tellheight=false
    )






    spaghetti_toggle = Toggle(mode_layout[2, 1], active=false)
    hexbin_toggle = Toggle(mode_layout[2, 2], active=true)
    surface_toggle = Toggle(mode_layout[2, 3], active=true)

    contours_observable = @lift([contours(eof_result.lon, eof_result.lat, all_correlation_data[$mode][member_id, :, :, $current_scope_index], $contour_levels) for member_id in 1:nmember])

    # contour_lines_observable = @lift(get_isocontour_lines($contours_observable...))

    seperate_matrices_observables = [@lift(all_correlation_data[$mode][i, :, :, $current_scope_index]) for i in 1:nmember]

    surf = surface!(
        axis,
        lonmin .. lonmax,
        latmin .. latmax,
        @lift(dropdims(mean(all_correlation_data[$mode][:, :, :, $current_scope_index], dims=1), dims=1));
        shading=shading,
        colormap=colormap,
        colorrange=limits,
        overdraw=false,
        transformation=(; translation=(0, 0, -1 * limits[2]))
    )

    connect!(surf.visible, surface_toggle.active)

    lines!(axis, GeoMakie.coastlines(); color=coastline_color, transformation=(; translation=(0, 0, 1000)))

    spaghetti_contours = [contour!(
        axis,
        lonmin .. lonmax,
        latmin .. latmax,
        matrix,
        levels=contour_levels,
        color=contour_colors[i],
        transformation=(; translation=(0, 0, 1100))
    ) for (i, matrix) in enumerate(seperate_matrices_observables)]

    for c in spaghetti_contours
        connect!(c.visible, spaghetti_toggle.active)
    end

    ensemble_vertices_obs = @lift([vcat(sample_along_line.(lines; dx=sample_distance)...) for lines in get_isocontour_lines($contours_observable...)])

    hexbin_plot = ensemblehexbin!(axis, ensemble_vertices_obs, cellsize=2.0, colormap=hexbin_colormap, threshold=1)

    Colorbar(layout[1, 3], hexbin_plot, label="Member Having Contour Line in Bin", vertical=true)


    picontrol_observable_correltaion = @lift(get_correlation_map(piControl.members[1].data[:, :, all_scopes[$current_scope_index]], eof_result.piControl[$current_scope_index].temporal_modes[:, $mode]))

    contour!(
        axis,
        lonmin .. lonmax,
        latmin .. latmax,
        picontrol_observable_correltaion,
        levels=contour_levels,
        color=:red,
        transformation=(; translation=(0, 0, 1100))
    )




    connect!(hexbin_plot.visible, hexbin_toggle.active)



    Colorbar(layout[1, 2], limits=limits, colormap=colormap, vertical=true, label=data_label)


end






function modes_stat_analysis(ensembles::Tuple{EOFEnsemble,String}...; resolution=(1080, 1920), fontsize=12, al_position=:rb, nmodes=nothing, stat_fun_type="variance", unit="")



    fig = Figure(size=resolution, fontsize=fontsize)


    stat_fun = lowercase(stat_fun_type) == "variance" ? var : std

    stat_label = lowercase(stat_fun_type) == "variance" ? "Variance" : "Standard Deviation"

    println("Calculates $stat_label")

    all_axis = []

    for (i, (eof_ensemble, data_name)) in enumerate(ensembles)

        scopes = eof_ensemble.scopes

        time_axis = eof_ensemble.time

        if isnothing(nmodes)
            nmodes = length(get_modes_variability(eof_ensemble.ensemble[get_member_id_string(11)][1]))
        end

        xmappings = Int[]

        ymappings = [Float64[] for _ in 1:nmodes]

        piControly = [Float64[] for _ in 1:nmodes]

        for scope_index in eachindex(scopes)

            for (_, eof_results) in eof_ensemble.ensemble

                results = stat_fun(eof_results[scope_index].temporal_modes, dims=1)

                push!(xmappings, scope_index)

                for mode in 1:nmodes
                    push!(ymappings[mode], results[mode])
                end

            end

            piControl_var = stat_fun(eof_ensemble.piControl[scope_index].temporal_modes, dims=1)

            for mode in 1:nmodes
                push!(piControly[mode], piControl_var[mode])
            end
        end

        axis = Axis(fig[i, 1], title="$(stat_label) of EOF coefficients of $data_name", ylabel=unit)

        push!(all_axis, axis)



        for mode in 1:nmodes
            boxplot!(axis, xmappings, ymappings[mode], label="Mode $mode")
        end
        for mode in 1:nmodes
            lines!(axis, eachindex(scopes), piControly[mode], color=:red, label="piControl")
        end



        seasonslice = 1:3:length(scopes)

        axis.xticks = (seasonslice, ["$(year(time_axis[scopes[i].start])) - $(year(time_axis[scopes[i].stop]))" for i in seasonslice])
        axis.xticklabelrotation = π / 4
        axis.xticklabelalign = (:right, :center)


        axislegend(axis, position=al_position, unique=true)
    end

    linkyaxes!(all_axis...)



    return fig
end
using EmpiricalOrthogonalFunctions
using NCDatasets
using NetCDF
using Dates
using BenchmarkTools
using Statistics
using PythonCall
using Distributed
using JLD2
using ProgressBars

include("utils.jl")


function get_eof_of_datachunk(datachunk; nmodes=nothing, center=true, alignment_field=nothing, varimax_rotation=false, scale_with_eigenvals=false)::EOFResult

    eof = EmpiricalOrthogonalFunction(datachunk; timedim=3, center=center)

    if isnothing(nmodes)
        nmodes = size(datachunk, 3)
    end

    if varimax_rotation
        orthorotation!(eof; n=nmodes)
    end

    temporalsignal = pcs(eof; n=nmodes)
    spatialsignal = eofs(eof; n=nmodes)

    if !isnothing(alignment_field)
        spatialsignal = align_with_field(spatialsignal, alignment_field)
    end

    if scale_with_eigenvals
        spatialsignal = spatialsignal .* sqrt.(eigenvalues(eof; n=nmodes))
    end

    modes_variability = eigenvalues(eof; n=nmodes) ./ sum(eof.eigenvals) * 100
    return EOFResult(reshape(spatialsignal, (size(datachunk)[1:2]..., nmodes)), temporalsignal, sqrt.(eigenvalues(eof; n=nmodes)), sum(eof.eigenvals), noscaling)
end

function pyeof_of_datachunk(datachunk, nmodes; weights=nothing, standard_permute=true, eof_type=:normal, eof_alignment_field=nothing, pcs_alignment_field=nothing, scale_mode=nothing, ddof=1.0)


    Eof = pyimport("eofs.standard").Eof
    np = pyimport("numpy")

    # eofs expects time to be the first dimension
    if standard_permute
        correct_shape_ds = permutedims(datachunk, (3, 1, 2))
    else
        correct_shape_ds = datachunk
    end

    if isnothing(scale_mode)
        scaling = 0
    elseif scale_mode == :divide
        scaling = 1
    elseif scale_mode == :multiply
        scaling = 2
    elseif scale_mode == :singularvals
        scaling = 0
    end


    if !isnothing(weights)
        solver = Eof(np.array(correct_shape_ds, dtype=np.float64), weights=np.array(weights, dtype=np.float64), ddof=ddof)
    else
        solver = Eof(np.array(correct_shape_ds, dtype=np.float64), ddof=ddof)
    end

    if eof_type == :normal
        eof_res = solver.eofs(neofs=nmodes, eofscaling=scaling)
    elseif eof_type == :covariance
        eof_res = solver.eofsAsCovariance(neofs=nmodes)
    elseif eof_type == :correlation
        eof_res = solver.eofsAsCorrelation(neofs=nmodes)
    else
        ArgumentError("Could not use eof type $(eof_type), use :normal, :correlation or :covariance instead")
    end

    modes_variability = solver.eigenvalues(neigs=nmodes) ./ solver.totalAnomalyVariance() * 100


    # we expect time to be the last dimension
    if !isnothing(eof_alignment_field)
        aligned_res = align_with_field(permutedims(pyconvert(Array{Float64,3}, eof_res), (2, 3, 1)), eof_alignment_field; mode_dim_index=2)
        spatialsignal = reshape(aligned_res, (size(datachunk)[1:2]..., nmodes))
    else
        spatialsignal = permutedims(pyconvert(Array{Float64,3}, eof_res), (2, 3, 1))
    end

    if !isnothing(pcs_alignment_field)
        aligned_pcs_res = align_with_field(pyconvert(Matrix{Float64}, solver.pcs(npcs=nmodes, pcscaling=scaling)), pcs_alignment_field; mode_dim_index=2)
        temporalsignal = reshape(aligned_pcs_res, (size(datachunk)[3], nmodes))
    else
        temporalsignal = pyconvert(Matrix{Float64}, solver.pcs(npcs=nmodes, pcscaling=scaling))
    end

    # revert the change from singular values to eigenvalues here https://github.com/ajdawson/eofs/blob/main/lib/eofs/standard.py
    if scale_mode == :singularvals
        eigenvals = pyconvert(Vector{Float64}, solver.eigenvalues(neigs=nmodes))
        singular_vals = sqrt.(eigenvals * (size(datachunk, 3) - ddof))
        spatialsignal = cat(
            [spatialsignal[:, :, mode] .* singular_vals[mode] * (1 / sqrt(size(datachunk, 3) - 1)) for mode in 1:nmodes]...,
            dims=3
        )
    end

    EOFResult(spatialsignal, temporalsignal, pyconvert(Vector{Float64}, solver.eigenvalues(neigs=nmodes)), pyconvert(Float64, solver.totalAnomalyVariance()), noscaling)

end


function calculate_eofs_of_tl_data(tl_data::TimelineData, chunking, nmodes; engine=:julia, reof=false, center=true, align_eof_with_mean=true, align_pcs_with_mean=false, weights=nothing, eof_type=:normal, scale_mode=nothing)::Dict{String,Vector{EOFResult}}

    function calculate_eof(chunk; eof_alignment_field=nothing, pcs_alignment_field=nothing)
        # if engine == :julia
        #     return get_eof_of_datachunk(chunk; nmodes=nmodes, center=center, alignment_field=alignment_field, varimax_rotation=reof, scale_with_eigenvals=scale_with_eigenvals)
        # elseif engine == :python
        #     return pyeof_of_datachunk(chunk, nmodes; weights=weights, standard_permute=true, eof_type=eof_type, eof_alignment_field=eof_alignment_field, pcs_alignment_field=pcs_alignment_field, scale_mode=scale_mode)
        # else
        #     ArgumentError("Could not recognize engine. Please use :python or :julia")
        # end
        return pyeof_of_datachunk(chunk, nmodes; weights=weights, standard_permute=true, eof_type=eof_type, eof_alignment_field=eof_alignment_field, pcs_alignment_field=pcs_alignment_field, scale_mode=scale_mode)
    end

    result = Dict()

    for (i, dataset) in enumerate(tl_data.datasets)
        eofs = Vector{EOFResult}(undef, length(chunking))  # Predefined array for EOFResult

        # Function to handle EOF calculation
        function handle_eof(idx, data, scopes, eofs)
            scope = scopes[idx]
            chunk = data[:, :, scope]

            eof_alignment_field = align_eof_with_mean ? reshape(mean(chunk, dims=3), 1, :) : nothing
            pcs_alignment_field = align_pcs_with_mean ? mean(chunk, dims=[1, 2]) : nothing
            # eof_alignment_field = align_eof_with_mean ? ones(prod(size(chunk)[1:2])) : nothing
            # pcs_alignment_field = align_pcs_with_mean ? ones(size(chunk)[3]) : nothing
            eof_result = calculate_eof(chunk; eof_alignment_field=eof_alignment_field, pcs_alignment_field=pcs_alignment_field)

            println("Handled scope $scope out of $(length(scopes)) on thread $(Threads.threadid())")

            # Direct assignment to the predefined array position
            eofs[idx] = eof_result
        end

        # Decide whether to use threading based on CONDITION
        if engine == :julia
            Threads.@threads for idx in eachindex(chunking)
                handle_eof(idx, dataset.data, chunking, eofs)
            end
        elseif engine == :python
            for idx in eachindex(chunking)
                handle_eof(idx, dataset.data, chunking, eofs)
            end
        else
            ArgumentError("COuld not recognize engine. Please use :python or :julia")
        end

        result[dataset.name] = eofs
    end

    return result

end

function calculate_eofs_of_ensemble(ensemble::EnsembleSimulation, chunking, nmodes; engine=:julia, reof=false, center=true, align_eof_with_mean=true, align_pcs_with_mean=false, weights=nothing, eof_type=:normal, scale_mode=nothing, saving_filepath=nothing)::Dict{String,Vector{EOFResult}}

    function calculate_eof(chunk; eof_alignment_field=nothing, pcs_alignment_field=nothing)
        # if engine == :julia
        #     return get_eof_of_datachunk(chunk; nmodes=nmodes, center=center, alignment_field=alignment_field, varimax_rotation=reof, scale_with_eigenvals=scale_with_eigenvals)
        # elseif engine == :python
        #     return pyeof_of_datachunk(chunk, nmodes; weights=weights, standard_permute=true, eof_type=eof_type, eof_alignment_field=eof_alignment_field, pcs_alignment_field=pcs_alignment_field, scale_with_eigenvals=scale_with_eigenvals)
        # else
        #     ArgumentError("Could not recognize engine. Please use :python or :julia")
        # end

        return pyeof_of_datachunk(chunk, nmodes; weights=weights, standard_permute=true, eof_type=eof_type, eof_alignment_field=eof_alignment_field, pcs_alignment_field=pcs_alignment_field, scale_mode=scale_mode)
    end

    result = Dict{String,Vector{EOFResult}}()

    function handle_eof(chunk)

        eof_alignment_field = align_eof_with_mean ? reshape(mean(chunk, dims=3), 1, :) : nothing
        pcs_alignment_field = align_pcs_with_mean ? mean(chunk, dims=[1, 2]) : nothing
        # eof_alignment_field = align_eof_with_mean ? ones(prod(size(chunk)[1:2])) : nothing
        # pcs_alignment_field = align_pcs_with_mean ? ones(size(chunk)[3]) : nothing
        eof_result = calculate_eof(chunk; eof_alignment_field=eof_alignment_field, pcs_alignment_field=pcs_alignment_field)

        # Direct assignment to the predefined array position
        return eof_result
    end

    for (i, member) in enumerate(ensemble.members)
        # Predefined array for EOFResult

        # Function to handle EOF calculation


        # Decide whether to use threading based on CONDITION
        @time "Time it took for eof calculation for member $(member.id)" begin
            if engine == :julia
                eofs = Vector{EOFResult}(undef, length(chunking))
                Threads.@threads for idx in eachindex(chunking)
                    eofs[idx] = handle_eof(chunk)
                end
            elseif engine == :python
                eofs = pmap([member.data[:, :, scope] for scope in chunking]) do chunk
                    handle_eof(chunk)
                end
            else
                ArgumentError("Could not recognize engine. Please use :python or :julia")
            end

            result[member.id] = eofs
        end

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
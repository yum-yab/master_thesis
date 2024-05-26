using LinearAlgebra
using Statistics
using EmpiricalOrthogonalFunctions


@enum EOFScaling begin
    noscaling = 1
    svalmult = 2
    svaldivide = 3
end

struct EOFResult
    spatial_modes::Array{Float64,3}
    temporal_modes::Array{Float64,2}
    singularvals::Array{Float64,1}
    sum_all_eigenvals::Float64
    scaling::EOFScaling
end


function scale_eof_result(eof_result::EOFResult; scale_mode::EOFScaling=svalmult)
    # now scale the eofs in different ways

    if eof_result.scaling != noscaling
        ArgumentError("You cannot scale an already scaled response!")
    end


    if scale_mode == svalmult
        new_eofs = eof_result.spatial_modes .* reshape(eof_result.singularvals, 1, 1, :)
        new_temp_signal = eof_result.temporal_modes .* reshape(eof_result.singularvals, 1, :)
    elseif scale_mode == svaldivide
        new_eofs = eof_result.spatial_modes ./ reshape(eof_result.singularvals, 1, 1, :)
        new_temp_signal = eof_result.temporal_modes ./ reshape(eof_result.singularvals, 1, :)
    end

    return EOFResult(new_eofs, new_temp_signal, eof_result.singularvals, eof_result.sum_all_eigenvals, scale_mode)
end

function get_modes_variability(eof_result::EOFResult)::Vector{Float64}
    return (eof_result.singularvals .^ 2) / eof_result.sum_all_eigenvals
end

function truncate_eof_result(eof_result::EOFResult, n::Int)
    return EOFResult(
        eof_result.spatial_modes[:, :, 1:n],
        eof_result.temporal_modes[:, 1:n],
        eof_result.singularvals[1:n],
        eof_result.sum_all_eigenvals,
        eof_result.scaling
    )
end

function reconstruct_data(eof_response::EOFResult, original_data::AbstractArray{<:Union{Missing,AbstractFloat},3}; timedim::Int=3, weights=nothing, center=true)


    # regenerate L, S, R depending on scaling
    S = Diagonal(eof_response.singularvals)

    # reshape spatial dims back
    # first two dims are geo_dims
    geo_shape = size(eof_response.spatial_modes)[1:2]

    reshaped_spatial = reshape(eof_response.spatial_modes, prod(geo_shape), :)

    if eof_response.scaling == noscaling
        L = eof_response.temporal_modes
        R = reshaped_spatial

    else
        broadcasting_svals = reshape(eof_response.singularvals, 1, :)
        if eof_response.scaling == svalmult
            L = eof_response.temporal_modes ./ broadcasting_svals
            R = reshaped_spatial ./ broadcasting_svals
        elseif eof_response.scaling == svaldivide
            L = eof_response.temporal_modes .* broadcasting_svals
            R = reshaped_spatial .* broadcasting_svals
        else
            ArgumentError("Unknown scaling: $(eof_response.scaling)")
        end
    end

    reconstructed = L * S * R'

    # reshape for timedim is the first, then permute the new dims 
    reconstructed = permutedims(reshape(reconstructed, :, geo_shape...), (2, 3, 1))


    if !isnothing(weights)
        weights_reshaped = reshape(weights, 1, :, 1)
        reconstructed = reconstructed ./ weights_reshaped
    end
    if center
        old_data_mean = mean(original_data, dims=timedim)

        # now reshape the reconstructed back to its original dimensions
        reconstructed .= reconstructed .+ old_data_mean
    end


    return reconstructed
end




function prepare_data_for_eof(data; weights=nothing, timedim=3, weightdim=2, norm_withsqrt_timedim=false)
    dims = size(data)

    time_shape = dims[timedim]

    geodims = [i for i in 1:length(dims) if i != timedim]
    geo_shape = dims[geodims]

    norm_factor = norm_withsqrt_timedim ? sqrt(time_shape - 1) : 1

    # apply weights over weightdim
    if !isnothing(weights)

        weights_reshaped = reshape(weights, 1, dims[weightdim], 1)
        data = data .* weights_reshaped ./ norm_factor
    end
    # reshape to timedim being the first
    if timedim != 1
        data = permutedims(data, (timedim, geodims...))
    end

    # flatten the data

    return reshape(data, (time_shape, prod(geo_shape)))
end

function align_with_field(field, alignment_field; mode_dim_index=2)

    dims = size(field)
    dimension = length(dims)

    if dimension == 3
        spat_dims = dims[1:2]
        mode_dim = dims[3]
        field = reshape(field, (prod(spat_dims), mode_dim))
    end

    function handle_slice(slice)
        scalar_product = sum(slice .* alignment_field)

        if scalar_product < 0
            return -1, slice * -1
        else
            return 1, slice
        end
    end

    result = similar(field)

    flips = ones(length(axes(field, mode_dim_index)))

    for modenum in axes(field, mode_dim_index)
        factor, sliceresult = handle_slice(field[:, modenum])
        result[:, modenum] = sliceresult
        flips[modenum] = factor
    end

    return result, flips
end


function eof(data; weights=nothing, timedim=3, weightdim=2, center=true, nmodes=-1, norm_withsqrt_timedim=false, align_eofs_with_mean=false)

    old_dims = size(data)

    geodims = [i for i in 1:length(old_dims) if i != timedim]

    data = prepare_data_for_eof(data; weights=weights, timedim=timedim, weightdim=weightdim, norm_withsqrt_timedim=norm_withsqrt_timedim)

    mean_over_time = mean(data, dims=1)

    time_shape = size(data, 1)

    if center
        data = data .- mean_over_time
    end

    if nmodes < 0
        mode_selector = 1:time_shape
    else
        mode_selector = 1:nmodes
    end
    # filter missing data

    left_singular_vectors, singular_vals, right_singular_vectors = svd(data)

    truncated_lsv = left_singular_vectors[:, mode_selector]
    truncated_svals = singular_vals[mode_selector]
    truncated_rsv = right_singular_vectors[:, mode_selector]

    # eigenvals
    eigenvals = singular_vals .^ 2

    eofs = truncated_rsv


    pc_factor = norm_withsqrt_timedim ? sqrt(time_shape - 1) : 1
    temporal_modes = truncated_lsv .* pc_factor



    if align_eofs_with_mean
        eofs, flip_factors = align_with_field(eofs, mean_over_time)
        temporal_modes = temporal_modes .* reshape(flip_factors, 1, :)
    end

    return EOFResult(reshape(eofs, (old_dims[geodims]..., :)), temporal_modes, truncated_svals, sum(eigenvals), noscaling)
end

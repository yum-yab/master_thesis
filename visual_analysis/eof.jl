using LinearAlgebra
using Statistics
using EmpiricalOrthogonalFunctions


@enum EOFScaling begin
    noscaling = 1
    tempmodescale = 2
    spatmodescale = 3
end

struct EOFResult
    spatial_modes::Array{Float64,3}
    temporal_modes::Array{Float64,2}
    singularvals::Array{Float64,1}
    sum_all_eigenvals::Float64
    scaling::EOFScaling
end


function scale_eof_result(eof_result::EOFResult; scale_mode::EOFScaling=spatmodescale)::EOFResult
    # now scale the eofs in different ways

    if eof_result.scaling != noscaling
        ArgumentError("You cannot scale an already scaled response!")
    end

    if scale_mode == spatmodescale 
        broadcast_svals = reshape(eof_result.singularvals, 1, 1, :)
        return EOFResult(eof_result.spatial_modes .* broadcast_svals, eof_result.temporal_modes, eof_result.singularvals, eof_result.sum_all_eigenvals, scale_mode)
    elseif scale_mode == tempmodescale
        broadcast_svals = reshape(eof_result.singularvals, 1, :)
        return EOFResult(eof_result.spatial_modes, eof_result.temporal_modes .* broadcast_svals, eof_result.singularvals, eof_result.sum_all_eigenvals, scale_mode)
    end
end


function scale_eof_result!(eof_result::EOFResult; scale_mode::EOFScaling=spatmodescale)
    # now scale the eofs in different ways

    if eof_result.scaling != noscaling
        ArgumentError("You cannot scale an already scaled response!")
    end

    if scale_mode == spatmodescale 
        broadcast_svals = reshape(eof_result.singularvals, 1, 1, :)
        eof_result.spatial_modes .*= broadcast_svals
    elseif scale_mode == tempmodescale
        broadcast_svals = reshape(eof_result.singularvals, 1, :)
        eof_result.temporal_modes .*= broadcast_svals
    end

    eof_result.scaling = scale_mode
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



function realign_modes!(eof_result::EOFResult, alignment_fields::Vector{Union{Nothing, Vector{<:Union{Missing,AbstractFloat}}}})
    
    for (i, align_field) in enumerate(alignment_fields)

        if isnothing(align_field)
            continue
        end

        factor = check_alignment_of_fields(eof_result.spatial_modes[:, :, i], align_field)

        eof_result.spatial_modes[:, :, i] *=  factor
        eof_result.temporal_modes[:, i] *=  factor
    end


end

function realign_modes(eof_result::EOFResult, alignment_fields)::EOFResult

    new_spatial_modes = similar(eof_result.spatial_modes)
    new_temporal_modes = similar(eof_result.temporal_modes)
    
    for (i, align_field) in enumerate(alignment_fields)

        if isnothing(align_field)
            continue
        end

        factor = check_alignment_of_fields(eof_result.spatial_modes[:, :, i], align_field)

        new_spatial_modes[:, :, i] = eof_result.spatial_modes[:, :, i] .* factor
        new_temporal_modes[:, i] = eof_result.temporal_modes[:, i] .* factor
    end

    return EOFResult(new_spatial_modes, new_temporal_modes, eof_result.singularvals, eof_result.sum_all_eigenvals, eof_result.scaling)
end

function reconstruct_data(eof_response::EOFResult; original_data_timemean::Union{Nothing,AbstractArray{<:Union{Missing,AbstractFloat},3}}=nothing)


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
        if eof_response.scaling == tempmodescale
            L = eof_response.temporal_modes ./ broadcasting_svals
            R = reshaped_spatial
        elseif eof_response.scaling == spatmodescale
            L = eof_response.temporal_modes
            R = reshaped_spatial ./ broadcasting_svals
        else
            ArgumentError("Unknown scaling: $(eof_response.scaling)")
        end
    end

    reconstructed = L * S * R'

    # reshape for timedim is the first, then permute the new dims 
    reconstructed = permutedims(reshape(reconstructed, :, geo_shape...), (2, 3, 1))

    if !isnothing(original_data_timemean)
        # now reshape the reconstructed back to its original dimensions
        reconstructed .= reconstructed .+ original_data_timemean
    end


    return reconstructed
end




function prepare_data_for_eof(data; meanfield=nothing, weights=nothing, timedim=3, weightdim=2, norm_withsqrt_timedim=false)
    dims = size(data)

    time_shape = dims[timedim]

    geodims = [i for i in 1:length(dims) if i != timedim]
    geo_shape = dims[geodims]

    # generate the needed factors
    norm_factor = norm_withsqrt_timedim ? 1 / sqrt(time_shape - 1) : 1
    data_adjusted = isnothing(meanfield) ? data : data .- meanfield
    weights_reshaped = isnothing(weights) ? ones(1, dims[weightdim], 1) : reshape(weights, 1, dims[weightdim], 1)

    # apply the needed factors
    data = norm_factor * weights_reshaped .* data_adjusted

    # reshape to timedim being the first
    if timedim != 1
        data = permutedims(data, (timedim, geodims...))
    end

    # flatten the data

    return reshape(data, (time_shape, prod(geo_shape)))
end

function check_alignment_of_fields(spatial_field, alignment_field)::Int

    scalar_product = dot(spatial_field, alignment_field)

    if scalar_product < 0
        return -1
    else
        return 1
    end
end

function align_all_modes_with_field(field, alignment_field; mode_dim_index=2)

    dims = size(field)
    dimension = length(dims)

    if dimension == 3
        spat_dims = dims[1:2]
        mode_dim = dims[3]
        field = reshape(field, (prod(spat_dims), mode_dim))
    end

    result = similar(field)

    flips = ones(length(axes(field, mode_dim_index)))

    for modenum in axes(field, mode_dim_index)
        factor = check_alignment_of_fields(field[:, modenum], alignment_field)
        result[:, modenum] = field[:, modenum] * factor
        flips[modenum] = factor
    end

    return result, flips
end

"""
    eof(data; weights=nothing, timedim=3, weightdim=2, center=true, nmodes=-1, norm_withsqrt_timedim=false, align_eofs_with_mean=false)

Compute the *nmodes* empirical orthogonal functions.  

If `y` is unspecified, compute the Bar index between all pairs of columns of `x`.

"""
function eof(data; weights=nothing, timedim=3, weightdim=2, center=true, nmodes=-1, norm_withsqrt_timedim=false, align_eofs_with_mean=false)

    old_dims = size(data)

    geodims = [i for i in 1:length(old_dims) if i != timedim]

    spat_dims = old_dims[geodims]

    other_dim = [i for i in 1:length(old_dims) if i != timedim && i != weightdim][1]

    if center || align_eofs_with_mean
        timemean = mean(data, dims=timedim)
    else
        timemean = nothing
    end

    if center
        mfield = timemean
    else
        mfield = nothing
    end

    data = prepare_data_for_eof(data; meanfield=mfield, weights=weights, timedim=timedim, weightdim=weightdim, norm_withsqrt_timedim=norm_withsqrt_timedim)

    time_shape = size(data, 1)

    if nmodes < 0
        mode_selector = Colon()
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
        weighted_timemean = timemean .* reshape(weights, 1, old_dims[weightdim], 1)
        # use the deviations from the spatial mean to compute the alignment, not the raw temporal mean
        eofs, flip_factors = align_all_modes_with_field(eofs, reshape(weighted_timemean .- mean(weighted_timemean), (prod(spat_dims), 1)))
        temporal_modes = temporal_modes .* reshape(flip_factors, 1, :)
    end

    if !isnothing(weights)
        eofs = eofs ./ repeat(weights, inner=(old_dims[other_dim], 1))
    end


    return EOFResult(reshape(eofs, (old_dims[geodims]..., :)), temporal_modes, truncated_svals, sum(eigenvals), noscaling)
end


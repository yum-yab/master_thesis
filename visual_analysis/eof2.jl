using LinearAlgebra
using Statistics


struct EOFResult
    spatial_modes::Array{Float64,3}
    temporal_modes::Array{Float64,2}
    modes_variability::Array{Float64,1}
end


function prepare_data_for_eof(data; weights=nothing, timedim=3, weightdim=2, norm_withsqrt_timedim=true)
    dims = size(data)

    time_shape = dims[timedim]

    geodims = [i for i in 1:length(dims) if i != timedim]
    geo_shape = dims[geodims]

    if norm_withsqrt_timedim

        data = data / sqrt(time_shape - 1)
    end

    # apply weights over weightdim
    if !isnothing(weights)

        weights_reshaped = reshape(weights, 1, dims[weightdim], 1)
        data = data .* weights_reshaped
    end

    # reshape to timedim being the first
    if timedim != 1
        data = permutedims(data, (timedim, geodims...))
    end

    # flatten the data

    return reshape(data, (time_shape, prod(geo_shape)))
end


function eof(data; weights=nothing, timedim=3, weightdim=2, center=true, scaling=nothing, nmodes=-1, norm_withsqrt_timedim=true)

    old_dims = size(data)

    geodims = [i for i in length(old_dims) if i != timedim]

    data = prepare_data_for_eof(data; weights=weights, timedim=timedim, weightdim=weightdim, norm_withsqrt_timedim=norm_withsqrt_timedim)

    time_shape = size(data, 1)

    geo_shape = size(data, 2)

    mean_over_time = mean(data, dims=1)

    if center
        data = data .- mean_over_time
    end

    if nmodes < 0
        mode_selector = 1:time_shape
    else
        mode_selector = 1:nmodes
    end
    # filter missing data

    isnotmissing_indices = findall(x -> !ismissing(x), mean_over_time[1, :])

    left_singular_vectors, singular_vals, right_singular_vectors = svd(data[:, isnotmissing_indices])

    truncated_lsv = left_singular_vectors[:, mode_selector]
    truncated_svals = singular_vals[mode_selector]
    truncated_rsv = right_singular_vectors[:, mode_selector]

    # eigenvals
    eigenvals = singular_vals .^ 2
    truncated_evs = eigenvals[mode_selector]

    # prepare svals to be broadcasted over modes
    broadcasting_svals = reshape(truncated_svals, 1, length(truncated_svals))

    # also reintroduce missing data

    eofs = Array{Union{Missing,Float64}}(missing, geo_shape, length(mode_selector))
    eofs[isnotmissing_indices, :] = truncated_rsv 


    # rescale the principal components with the singular value
    pc_factor = norm_withsqrt_timedim ? sqrt(time_shape - 1) : 1
    principal_components = truncated_lsv .* broadcasting_svals * pc_factor

    # now scale the eofs in different ways
    if !isnothing(scaling)

        if scaling == :singularvals

            eofs = eofs .* broadcasting_svals
        elseif scaling == :eigenvalsmult
            sqrt_ev = sqrt.(truncated_evs)
            eofs = eofs .* reshape(sqrt_ev, 1, length(sqrt_ev))
        elseif scaling == :eigenvalsdivide
            sqrt_ev = sqrt.(truncated_evs)
            eofs = eofs ./ reshape(sqrt_ev, 1, length(sqrt_ev))
        end

    end

    println(size(reshape(eofs, (old_dims[geodims]..., length(mode_selector)))))

    return EOFResult(reshape(eofs, (old_dims[geodims]..., length(mode_selector))), principal_components, truncated_evs / sum(eigenvals) * 100)
end
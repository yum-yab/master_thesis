include("generate_ivt_fields.jl")

using .preprocessing

using HDF5

function old_normal_generation(id_to_file_mapping, geo_bounds)::Nothing
    _ = generate_ivt_field(id_to_file_mapping, geo_bounds)
    return
end

function parallel_reading_at_start(id_to_file_mapping, geo_bounds)::Nothing

    ids = ["hus", "ua", "va", "ps"]

    acc = [[geo_bounds.lon_indices[1]:geo_bounds.lon_indices[1], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[1], Colon(), Colon()] for _ in 1:3]

    push!(acc, [geo_bounds.lon_indices[1]:geo_bounds.lon_indices[1], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[1], Colon()])

    data_lookup = Dict()

    lk = Threads.ReentrantLock()

    Threads.@threads for i in eachindex(ids)
        id = ids[i]
        path = id_to_file_mapping[id]
        h5open(path, "r") do h5ds
            data = h5ds[id][acc[i]...]
            lock(lk) do
                data_lookup[id] = data
            end
        end
    end

    ap = h5read(id_to_file_mapping["hus"], "ap")
    b = h5read(id_to_file_mapping["hus"], "b")

    result_data = zeros(Float32, lon_length, lat_length, time_length)


    Threads.@threads for time in 1:time_length

        for lon in 1:lon_length
            for lat in 1:lat_length
                ps = data_lookup["ps"][lon, lat, time]
                pressure_levels = ap + b * ps

                result_data[lon, lat, time] = IVT.ivt_of_column_vectors(
                    ps,
                    pressure_levels,
                    data_lookup["hus"][lon, lat, :, time],
                    data_lookup["ua"][lon, lat, :, time],
                    data_lookup["va"][lon, lat, :, time]
                )
            end
        end
    end

    return
end


function parallel_reading_each_timestep(id_to_file_mapping, geo_bounds)::Nothing
    file_dict = Dict([id => h5open(path, "r") for (id, path) in id_to_file_mapping])


    ap = file_dict["hus"]["ap"][:]
    b = file_dict["hus"]["b"][:]

    result_data = zeros(Float32, lon_length, lat_length, time_length)



    Threads.@threads for time in 1:$time_length
        hus_data = file_dict["hus"]["hus"][geo_bounds.lon_indices[1]:geo_bounds.lon_indices[1], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[1], :, time]
        ua_data = file_dict["ua"]["ua"][geo_bounds.lon_indices[1]:geo_bounds.lon_indices[1], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[1], :, time]
        va_data = file_dict["va"]["va"][geo_bounds.lon_indices[1]:geo_bounds.lon_indices[1], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[1], :, time]
        ps_data = file_dict["ps"]["ps"][geo_bounds.lon_indices[1]:geo_bounds.lon_indices[1], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[1], time]

        for lon in 1:lon_length
            for lat in 1:lat_length
                ps = ps_data[lon, lat]
                pressure_levels = ap + b * ps

                result_data[lon, lat, time] = IVT.ivt_of_column_vectors(ps, pressure_levels, hus_data[lon, lat, :], ua_data[lon, lat, :], va_data[lon, lat, :])
            end
        end
    end
    for id in ids
        close(file_dict[id])
    end
    return
end

function main()
    
    geo_bnds = 

    println("Time it took for parallel reading at start: ")


end
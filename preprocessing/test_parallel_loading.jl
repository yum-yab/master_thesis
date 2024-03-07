include("generate_ivt_fields.jl")

using .preprocessing

using HDF5

function old_normal_generation(id_to_file_mapping, geo_bounds)::Nothing
    _ = generate_ivt_field(id_to_file_mapping, geo_bounds)
    return
end

function parallel_reading_at_start(id_to_file_mapping, geo_bounds)::Nothing

    ids = ["hus", "ua", "va", "ps"]

    acc = [[geo_bounds.lon_indices[1]:geo_bounds.lon_indices[2], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[2], Colon(), Colon()] for _ in 1:3]

    push!(acc, [geo_bounds.lon_indices[1]:geo_bounds.lon_indices[2], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[2], Colon()])

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

    dimfile = h5open(id_to_file_mapping["hus"])
    ap = dimfile["ap"][:]
    b = dimfile["ap"][:]    
    close(dimfile)
    
    dims = size(data_lookup["hus"])
    time_length = dims[4]

    lon_length = dims[1]

    lat_length = dims[2]
    result_data = zeros(Float32, lon_length, lat_length, time_length)


    Threads.@threads for time in 1:time_length

        for lon in 1:lon_length
            for lat in 1:lat_length
                ps = data_lookup["ps"][lon, lat, time]
                pressure_levels = ap + b * ps

                result_data[lon, lat, time] = preprocessing.IVT.ivt_of_column_vectors(
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
    
    time_length = size(file_dict["hus"]["hus"], 4)

    lon_length = length(geo_bounds.remaining_lon_values)
    lat_length = length(geo_bounds.remaining_lat_values)

    ap = file_dict["hus"]["ap"][:]
    b = file_dict["hus"]["b"][:]

    result_data = zeros(Float32, lon_length, lat_length, time_length)



    Threads.@threads for time in 1:time_length
        hus_data = file_dict["hus"]["hus"][geo_bounds.lon_indices[1]:geo_bounds.lon_indices[2], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[2], :, time]
        ua_data = file_dict["ua"]["ua"][geo_bounds.lon_indices[1]:geo_bounds.lon_indices[2], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[2], :, time]
        va_data = file_dict["va"]["va"][geo_bounds.lon_indices[1]:geo_bounds.lon_indices[2], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[2], :, time]
        ps_data = file_dict["ps"]["ps"][geo_bounds.lon_indices[1]:geo_bounds.lon_indices[2], geo_bounds.lat_indices[1]:geo_bounds.lat_indices[2], time]

        for lon in 1:lon_length
            for lat in 1:lat_length
                ps = ps_data[lon, lat]
                pressure_levels = ap + b * ps

                result_data[lon, lat, time] = preprocessing.IVT.ivt_of_column_vectors(ps, pressure_levels, hus_data[lon, lat, :], ua_data[lon, lat, :], va_data[lon, lat, :])
            end
        end
    end
    for (_, fileid) in file_dict
        close(fileid)
    end
    return
end

function print_and_test(fun, name, args...)
  
  println("Time it took for $name: ")
  @time fun(args...)
end

function main()
  println("Available Threads: $(Threads.nthreads())")    
  lon_bnds = (0,130)
  lat_bnds = (20, 80)
  
  id_to_file_mapping = Dict(
    "hus" => "/work/ik1017/CMIP6/data/CMIP6/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp585/r1i1p1f1/6hrLev/hus/gn/v20190710/hus_6hrLev_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_201501010600-203501010000.nc",
    "ua" => "/work/ik1017/CMIP6/data/CMIP6/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp585/r1i1p1f1/6hrLev/ua/gn/v20190710/ua_6hrLev_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_201501010600-203501010000.nc",
    "va" => "/work/ik1017/CMIP6/data/CMIP6/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp585/r1i1p1f1/6hrLev/va/gn/v20190710/va_6hrLev_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_201501010600-203501010000.nc",
    "ps" => "/work/ik1017/CMIP6/data/CMIP6/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/ssp585/r1i1p1f1/6hrLev/ps/gn/v20190710/ps_6hrLev_MPI-ESM1-2-LR_ssp585_r1i1p1f1_gn_201501010600-203501010000.nc",
  )

  geo_bounds = preprocessing.DataLoading.GeographicBounds(lon_bnds, lat_bnds, id_to_file_mapping["hus"])

  # println("Running the old way once to compile it:")
  # old_normal_generation(id_to_file_mapping, geo_bounds)
  
  print_and_test(parallel_reading_each_timestep, "parallel reading each timestep", id_to_file_mapping, geo_bounds) 
  print_and_test(parallel_reading_at_start, "parallel reading at start", id_to_file_mapping, geo_bounds) 
  print_and_test(old_normal_generation, "old normal generation", id_to_file_mapping, geo_bounds) 

end

main()

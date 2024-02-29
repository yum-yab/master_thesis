
module preprocessing

using NCDatasets
using DataStructures
using NetCDF

include("IVT.jl")
using .IVT

include("data_loading.jl")
using .DataLoading

export generate_ivt_field, GeographicBounds, load_data_in_geo_bounds




function generate_ivt_field(hus_field_path::String, ua_field_path::String, va_field_path::String, ps_field_path::String, output_path::String, lon_bounds::Tuple{<:Real,<:Real}, lat_bounds::Tuple{<:Real,<:Real})::Nothing


  
  geo_bnds = GeographicBounds(lon_bounds, lat_bounds, hus_field_path)


  println("Time used for loading the data: ")
  @time begin
    hus_data = load_variable_data_in_bounds(Float32, hus_field_path, "hus", geo_bnds, :, :)
    ua_data = load_variable_data_in_bounds(Float32, ua_field_path, "ua", geo_bnds, :, :)
    va_data = load_variable_data_in_bounds(Float32, va_field_path, "va", geo_bnds, :, :)
    ps_data = load_variable_data_in_bounds(Float32, ps_field_path, "ps", geo_bnds, :)

    lon_size = size(hus_data, 1)
    lat_size = size(hus_data, 2)

    # these variables are used for calculation of pressure levels at each specific lat, lon, time coordinate: p = ap + b * ps
    # using netcdf lib here since it will be quite fast
    ap = ncread(hus_field_path, "ap")
    b = ncread(hus_field_path, "b")
    time_size = size(hus_data, 4)

  end

  result_data = zeros(Float32, lon_size, lat_size, time_size)
  println("Time used for calculating the IVT field: ")
  @time begin
    Threads.@threads for time in 1:time_size
      for lat in 1:lat_size
        for lon in 1:lon_size

          ps = ps_data[lon, lat, time]
          pressure_levels = ap + b * ps

          result_data[lon, lat, time] = IVT.ivt_of_column_vectors(ps, pressure_levels, hus_data[lon, lat, :, time], ua_data[lon, lat, :, time], va_data[lon, lat, :, time])

        end
      end
    end

  end


  dataset = NCDataset(hus_field_path)

  # copy the data to a new Ordered Dict since a direct pass would only pass the reference to the old atrrib dict -> ERROR!
  new_attrib = OrderedDict()

  for (k, v) in dataset.attrib
    new_attrib[k] = v
  end


  # then rewrite some fields
  new_attrib["history"] = "Used an IVT julia script found here: https://github.com/yum-yab/master_thesis/tree/main/preprocessing" * new_attrib["history"]
  new_attrib["title"] = "IVT field of europe and the north-east Atlantic"

  NCDataset(output_path, "c", attrib=new_attrib) do ds

    defVar(ds, "lon", geo_bnds.remaining_lon_values, ("lon",), attrib=OrderedDict(
      "units" => "degrees east",
      "standard_name" => "longitude"
    ))
    defVar(ds, "lat", geo_bnds.remaining_lat_values, ("lat",), attrib=OrderedDict(
      "units" => "degrees north",
      "standard_name" => "latitude"
    ))
    defVar(ds, "time", dataset[:time][:], ("time",), attrib=OrderedDict(
      "units" => "days since 1850-1-1 00:00:00",
      "standard_name" => "time"
    ))
    defVar(ds, "ivt", result_data, ("lon", "lat", "time"), attrib=OrderedDict(
      "units" => "kg/ms",
      "comments" => "This field is NOT normalized over time or temperature or whatever"
    ))


  end

  close(dataset)

  return
end

end

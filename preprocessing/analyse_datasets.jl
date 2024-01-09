using NetCDF


BASE_PATH = ""

function info_of_dataset(ssp_scenario:: String, member_id:: int, time_resolution:: String, variable:: String, only_first_file:: Bool)
  
  dataset_path = "$BASE_PATH/$ssp_scenario/r$(member_id)i1p1f1/"
end

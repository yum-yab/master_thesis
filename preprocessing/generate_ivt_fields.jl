using NetCDF



function prepare_merged_file(input_files::Vector{String}, output_file::String)::String

  PRESSURE_LEVELS=round.(Int, LinRange(100000, 0, 47))
  
  LON_BOUNDS=(270, 40)
  LAT_BOUNDS=(20, 80)

  cdo_cmd = Cmd(["cdo", "-ml2pl,$(join(string.(PRESSURE_LEVELS), ","))", "-sellonlatbox,$(LON_BOUNDS[1]),$(LON_BOUNDS[2]),$(LAT_BOUNDS[1]),$(LAT_BOUNDS[2])", "-selvar,hus,ua,va", "-merge", input_files..., output_file])
  
  return run(cdo_cmd, wait=true)
end

function nc_files_in_directory(directory_path::String)::Vector{String}
  
  return filter(path -> endswith(path, ".nc"), readdir(directory_path, join = true))
end


function find_common_versions(member_path::String, field_ids::Vector{String}, time_res_id::String)::Union{String, Nothing}

  available_versions = [Set(readdir(joinpath(member_path, time_res_id, fid, "gn"))) for fid in field_ids]

  common_versions = Vector(intersect(available_versions...))

  if isempty(common_versions)
    return nothing 
  else
    # get the latest common version
    return sort(common_versions, rev = true)[1]
  end
end



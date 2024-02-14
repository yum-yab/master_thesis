import os
from typing import List, Optional, Tuple
import cdo
import math



CDO = cdo.Cdo()

def ncfiles_in_dir(directory: str) -> List[str]:
    
    result = [os.path.join(directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.endswith(".nc")] 
    result.sort()
    return result

def get_timestamp_from_name(filename: str) -> str:

    return filename.split("_")[-1].split(".")[0]


def generate_desired_pressure_levels(stepwidth: int, lims: Tuple[int, int] = (10000, 0)) -> List[int]:
    
    return [i for i in range(lims[0], lims[1], -1 * stepwidth)] + [lims[1]]



def generate_intermediate_pressure_levels(pressure_lvl_list: List[int]) -> List[float]:

    # assert list is sorted descending
    assert all(pressure_lvl_list[i] >= pressure_lvl_list[i+1] for i in range(len(pressure_lvl_list)-1))
    
    return [(pressure_lvl_list[i] + pressure_lvl_list[i-1])/2 for i in range(1,len(pressure_lvl_list)-1)]




def variable_merge_list_check(files_iterator):

    for filelist in files_iterator:
        timestamp = get_timestamp_from_name(filelist[0])
        for f in filelist:
            assert get_timestamp_from_name(f) == timestamp

def find_common_version(member_path: str, field_ids: List[str], time_res_id: str, ) -> Optional[str]:
    
    available_versions = [set(os.listdir(os.path.join(member_path, time_res_id, fid, "gn"))) for fid in field_ids]
    common_versions = list(set.intersection(*available_versions))
    
    if common_versions:
        common_versions.sort(reverse=True)
        return common_versions[0]
    else:
        return None


def run_cdo_ivt_field_gen(sourcefiles: List[str], target_file: str, pressure_levels: List[int]):

    p_lvl_string = ",".join([str(p) for p in pressure_levels])

    intermediate_p_lvl_string = ",".join([str(p) for p in generate_intermediate_pressure_levels(pressure_levels)])
    
    higher_vlv_indexes = ",".join([str(i) for i in list(range(1,len(pressure_levels)-1))])
    lower_vlv_indexes = ",".join([str(i) for i in list(range(0,len(pressure_levels)-2))])

    f"-intlevel,{intermediate_p_lvl_string} -ml2pl,{p_lvl_string} -select,ua,va,hus,ps {' '.join(sourcefiles)} {target_file}"





def generate_ivt_fields(data_base_path: str, scenario_id: str, target_base_path: str, time_res_id: str = "6hrLev", dry_run: bool = True, debug_cdo: bool = True, cdo_nthreads: Optional[int] = None):
    
    scenario_path = os.path.join(data_base_path, scenario_id)

    field_ids = ["hus", "ua", "va"]
    
    # members as (id, path) tuples
    members = [(d, os.path.join(scenario_path, d)) for d in os.listdir(scenario_path) if os.path.isdir(os.path.join(scenario_path, d))]
    
    cdo_inst = cdo.Cdo()

    cdo_default_options = ["-T"] if dry_run else []

    if cdo_nthreads:
        cdo_default_options.append(f"-P {str(cdo_nthreads)}")

    cdo_inst.debug = debug_cdo

    for member_id, member_path in members:

        # find available variable for member, if many are present take the latest run 
         
        version_id = find_common_version(member_path, field_ids, time_res_id)

        if version_id is None:
            print(f"No common version found for member {member_id} in scenario {scenario_id}!")
            continue
        else:
            print(f"Found common version for member {member_id} in scenario {scenario_id}: {version_id}")


        variable_paths = [os.path.join(member_path, time_res_id, variable, "gn", version_id) for variable in field_ids]
        # zips together each timestep from each variable for merging
        files_to_merge = list(zip(ncfiles_in_dir(variable_paths[0]), ncfiles_in_dir(variable_paths[1]), ncfiles_in_dir(variable_paths[2])))
        
        print(f"Files to merge {files_to_merge}")
        variable_merge_list_check(files_to_merge)

        for filelist in files_to_merge:
            
            timestamp = get_timestamp_from_name(filelist[0])
            
            target_path = os.path.join(target_base_path, scenario_id, member_id)

            os.makedirs(target_path, exist_ok=True)

            target_file = os.path.join(target_path, f"merged-fields_ua-va-hus_{scenario_id}_{member_id}_{timestamp}.nc")
            
            

            cdo_inst.merge(input = " ".join(filelist), output = target_file, options = " ".join(cdo_default_options))




def run_all_scps(data_base_path: str, target_base_path: str, dry_run: bool = True, debug_cdo: bool = True):

    ssp_ids = [d for d in os.listdir(data_base_path) if os.path.isdir(os.path.join(data_base_path, d))]

    for ssp in ssp_ids:

        generate_ivt_fields(data_base_path, ssp, target_base_path, dry_run=dry_run, debug_cdo=debug_cdo)


def main():

    dataset_dir = "/pool/data/CMIP6/data/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/"

    target_dir = "/scratch/b/b382641/merged_mpige_cmip6_data/"

    run_all_scps(dataset_dir, target_dir, dry_run=False)

if __name__ == "__main__":
    main()

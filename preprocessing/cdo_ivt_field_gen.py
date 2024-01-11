import subprocess
import os
from typing import List, Optional

def ncfiles_in_dir(directory: str) -> List[str]:
    
    result = [os.path.join(directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.endswith(".nc")] 
    result.sort()
    return result

def get_timestamp_from_name(filename: str) -> str:

    return filename.split("_")[-1].split(".")[0]


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


def merge_relevant_datasets(data_base_path: str, scenario_id: str, target_base_path: str, time_res_id: str = "6hrLev", dry_run: bool = True):
    
    scenario_path = os.path.join(data_base_path, scenario_id)

    field_ids = ["hus", "ua", "va"]
    
    # members as (id, path) tuples
    members = [(d, os.path.join(scenario_path, d)) for d in os.listdir(scenario_path) if os.path.isdir(os.path.join(scenario_path, d))]
    
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
        files_to_merge = zip(ncfiles_in_dir(variable_paths[0]), ncfiles_in_dir(variable_paths[1]), ncfiles_in_dir(variable_paths[2]))
        
        print(f"Files to merge {len(list(files_to_merge))}")
        variable_merge_list_check(files_to_merge)

        for filelist in files_to_merge:
            
            filepaths_str = " ".join(filelist)

            timestamp = get_timestamp_from_name(filelist[0])
            
            target_path = os.path.join(target_base_path, scenario_id, member_id)

            os.makedirs(target_path, exist_ok=True)

            target_file = os.path.join(target_path, f"merged-fields_ua-va-hus_{scenario_id}_{member_id}_{timestamp}.nc")

            command = ["cdo", "-A", "merge", f"[ {filepaths_str} ]", target_file] if dry_run else ["cdo", "merge", f"[ {filepaths_str} ]", target_file]
            # now call the merge command for one timestep
            print(f"CDO command: {command}")
            com_process = subprocess.run(command, encoding="utf-8", stderr=subprocess.STDOUT)
            
            print(com_process.stdout)

            com_process.check_returncode()


def run_all_scps(data_base_path: str, target_base_path: str, dry_run: bool = True):

    ssp_ids = [d for d in os.listdir(data_base_path) if os.path.isdir(os.path.join(data_base_path, d))]

    for ssp in ssp_ids:

        merge_relevant_datasets(data_base_path, ssp, target_base_path, dry_run=dry_run)


def main():

    dataset_dir = "/pool/data/CMIP6/data/ScenarioMIP/MPI-M/MPI-ESM1-2-LR/"

    target_dir = "/scratch/b/b382641/merged_mpige_cmip6_data/"

    run_all_scps(dataset_dir, target_dir, dry_run=True)

if __name__ == "__main__":
    main()

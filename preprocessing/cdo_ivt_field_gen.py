import subprocess
import os
from typing import Tuple 

DATA_BASE_PATH=""

def ncfiles_in_dir(directory: str):
    
    result = [f for f in os.listdir(directory) if os.path.isfile(f) and f.endswith(".nc")] 
    result.sort()
    return result

def get_timestamp_from_name(filename: str):

    return filename.split("_")[-1].split(".")[0]


def variable_merge_list_check(files_iterator):

    for filelist in files_iterator:
        timestamp = get_timestamp_from_name(filelist[0])
        for f in filelist:
            assert get_timestamp_from_name(f) == timestamp


def merge_relevant_datasets(scenario_id: str, target_base_path: str, time_res_id: str = "6hrLev", version_id: str = "v20190710"):
    
    scenario_path = os.path.join(DATA_BASE_PATH, scenario_id)
    
    # members as (id, path) tuples
    members = [(d, os.path.join(scenario_path, d)) for d in os.listdir(scenario_path) if os.path.isdir(os.path.join(scenario_path, d))]
    
    for member_id, member_path in members:

        variable_paths = [os.path.join(member_path, time_res_id, variable, "gn", version_id) for variable in ["hus", "ua", "va"]]
        
        # zips together each timestep from each variable for merging
        files_to_merge = zip(ncfiles_in_dir(variable_paths[0]), ncfiles_in_dir(variable_paths[1]), ncfiles_in_dir(variable_paths[2]))
        
        variable_merge_list_check(files_to_merge)

        for filelist in files_to_merge:

            filepaths_str = " ".join(filelist)

            timestamp = get_timestamp_from_name(filelist[0])
            
            target_path = os.path.join(target_base_path, scenario_id, member_id)

            os.makedirs(target_path, exist_ok=True)

            target_file = os.path.join(target_path, f"merged-fields_ua-va-hus_{scenario_id}_{member_id}_{timestamp}.nc")

            command = ["cdo", "merge", f"[ {filepaths_str} ]", target_file]
            # now call the merge command for one timestep 
            pass



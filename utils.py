import numpy as np

##### important data sets - support ####
bp_1 = np.array([1, 1, 1, 3, 3, 3, 5, 5, 5, 7, 7, 7])
bp_2 = np.array([0, 0, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7])
bp_3 = np.array([1, 1, 1, 2, 2, 3, 5, 5, 5, 7, 7, 7])
bp_4 = np.array([1, 1, 1, 3, 3, 3, 4, 4, 5, 7, 7, 7])
bp_5 = np.array([0, 0, 1, 2, 2, 3, 5, 5, 5, 6, 6, 7])
bp_6 = np.array([0, 0, 1, 2, 2, 3, 4, 4, 5, 7, 7, 7])
bp_7 = np.array([0, 0, 1, 3, 3, 3, 4, 4, 5, 7, 7, 7])
bp_8 = np.array([1, 1, 1, 3, 3, 3, 4, 4, 5, 6, 6, 7])
bp_9 = np.array([1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 7, 7])

bp_13 = np.array([1, 1, 1, 3, 3, 3, 5, 5, 5, 6, 6, 7])
bp_14 = np.array([0, 0, 1, 3, 3, 3, 5, 5, 5, 7, 7, 7])
bp_15 = np.array([1, 1, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7])
bp_16 = np.array([0, 0, 1, 3, 3, 3, 4, 4, 5, 6, 6, 7])
bp_17 = np.array([1, 1, 1, 2, 2, 3, 5, 5, 5, 6, 6, 7])
bp_18 = np.array([0, 0, 1, 2, 2, 3, 5, 5, 5, 7, 7, 7])
bp_19 = np.array([0, 0, 1, 3, 3, 3, 5, 5, 5, 6, 6, 7])

bp_dict = {
    'bp_1': bp_1,
    'bp_2': bp_2,
    'bp_3': bp_3,
    'bp_4': bp_4,
    'bp_5': bp_5,
    'bp_6': bp_6,
    'bp_7': bp_7,
    'bp_8': bp_8,
    'bp_9': bp_9,
    'bp_13': bp_13,
    'bp_14': bp_14,
    'bp_15': bp_15,
    'bp_16': bp_16,
    'bp_17': bp_17,
    'bp_18': bp_18,
    'bp_19': bp_19
}


save_final_results_filename = "run-1-results.txt"
save_all_raw_filename = ""

##### functions - support ####
# save ALL raw results
def save_all_raw(data_list):
    # this will save every bond_permute and output value (Decimal level precision) to a text file
    with open(save_all_raw_filename, 'w') as f:
        for bond_permute, result in data_list:
            f.write(f"{bond_permute} {result}\n")

def save_final_results(bp_label, start_time_str, end_time_str, runtime_minutes, cores, total_sum):
    with open(save_final_results_filename, 'a') as file:
        file.write(f"bp_label: {bp_label}\n")
        file.write(f"start time: {start_time_str}\n")
        file.write(f"end time: {end_time_str}\n")
        file.write(f"run time: {runtime_minutes:.2f} minutes\n")
        file.write(f"cores: {cores}\n")
        file.write(f"total_sum: {total_sum}\n")
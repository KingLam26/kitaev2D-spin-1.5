##### import #####
import numpy as np
import itertools


##### parameters: model #####
spin_S = 1.5
bonds_ext = np.array([0,2,4,6])
bonds_int = np.array([1,3,5,7])
y_bonds = [0,3,4,7]
bond_site_dict = {0: (1,0), 1: (1,9), 2: (3,2), 3: (3,9),
                  4: (5,4), 5: (5,8), 6: (7,6), 7: (7,8)}

y_bond_sites = [site for y_bond in y_bonds for site in bond_site_dict[y_bond]]

# ini_spin_config = [spin_S] * 12 # also use tuple..
# in parallel processing, avoid global mutable variables

gs_energy = spin_S * spin_S * 8


##### parameters: parallel processing #####
num_lockable_lists = 10
queue_maxsize = 10000
num_permutes_print = 25000
compute_switch = True
Trie_true = True

##### compute model parameters #####
def gen_A_combis_dict(combi_sizes):
    """
    Generate a dictionary of valid combinatorial arrays for given sizes.

    For each size in the list of `combi_sizes`, this function creates a base array 
    with an equal number of +1s and -1s (with possibly one extra +1 or -1 to balance 
    the array length). It then generates permutations of the middle part of this 
    base array and filters them based on a rolling sum condition. The result is a 
    dictionary where the key is the size of the combinatorial array and the value 
    is a NumPy array of valid permutations.

    Parameters:
    combi_sizes (list of int): List of sizes for which combinatorial arrays are to be generated.

    Returns:
    dict: A dictionary where each key is a size from `combi_sizes` and the corresponding 
          value is a NumPy array of a physical combinatorial arrays of that size.
    
    Example:
    >>> gen_A_combis_dict([4, 6])
    4 [[-1  1 -1  1]
        [-1  1  1 -1]
        [ 1 -1  1 -1]
        [ 1 -1 -1  1]
        [-1 -1  1  1]
        [ 1  1 -1 -1]]
    6 [[-1  1 -1  1 -1  1]
        [-1  1  1 -1 -1  1]
        [-1 -1  1  1 -1  1]
        [ 1 -1 -1  1  1 -1]
        ...
        [ 1 -1  1 -1  1 -1]]
    """


    A_combis_dict = {}
    
    for size in combi_sizes:
        # Create  base array with n/2 +1s and n/2 -1s
        half_n = size // 2
        base_array = [-1] + [1] * (half_n - 1) + [-1] * (half_n - 1) + [1]

        # Generate combis for the middle part
        middle_part = base_array[1:-1]
        combis = set(itertools.permutations(middle_part))

        # Filter permutations while generating final permutations
        filtered_combis = []
        for combi in combis:
            full_combi = list((-1,) + combi + (1,))
            roll_sum = np.cumsum(full_combi)
            if np.all(roll_sum <= 0):
                filtered_combis.append(full_combi)

        A_combis_dict[size] = np.array(filtered_combis)

    # print test
    """
    for key, value in A_combis_dict.items():
        print(key, value)
    """
    return A_combis_dict

A_combis_dict = gen_A_combis_dict(combi_sizes=[2,4,6])


##### print input parameters #####
print(f"num_lockable_lists: {num_lockable_lists}")
print(f"queue_maxsize: {queue_maxsize}")
print(f"num_permutes_print: {num_permutes_print}")
print(f"compute_switch: {str(compute_switch)}")

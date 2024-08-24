##### import #####
from inputs import *
import numpy as np
from collections import Counter
import math, time
from itertools import product

##### Trie class #####
class TrieNode:
    def __init__(self, spin_config, cum_energy_factor, sum_spin_factor):
        self.children = {}      # Dictionary where keys are labels (e.g., +1, -1)
        self.result = (spin_config, cum_energy_factor, sum_spin_factor)      # Store computation results for the current prefix


class Trie:
    def __init__(self, ini_spin_config, ini_energy_factor, ini_spin_factor):
        self.root = TrieNode(ini_spin_config, ini_energy_factor, ini_spin_factor)

    def insert_and_compute(self, A_site_combi, Z_sites, A_sites, combi_propagate):
        node = self.root

        for step in range(12):
            label = A_site_combi[step]  # Use the label (e.g., +1, -1) as a key

            if label not in node.children:
                # extract parameters from latest node
                spin_config, cum_energy_factor, cum_spin_factor = node.result   #tuple
                spin_config, local_energy_factor, local_spin_factor = combi_propagate(spin_config, step, Z_sites, A_sites, A_site_combi)

                if step != 11:
                    cum_energy_factor*=local_energy_factor
                cum_spin_factor*=local_spin_factor

                node.children[label] = TrieNode(spin_config, cum_energy_factor, cum_spin_factor)
            node = node.children[label] # move to child node

        cum_spin_factor = math.sqrt(cum_spin_factor)
        return cum_energy_factor, cum_spin_factor


##### functions #####
# generate permutes #
def next_unique_permute(array):

    """
    Generate the next lexicographical permutation of the input array in-place.

    This function modifies the input array to produce the next permutation in lexicographical order.
    If the array is already in the highest possible order, the function returns `False` and does not modify the array.
    Otherwise, it returns `True` after rearranging the array into the next permutation.

    Steps:
    1. Find the largest index `k` such that `array[k] < array[k + 1]`. If no such index exists, the array is the last permutation.
    2. Find the largest index `l` greater than `k` such that `array[k] < array[l]`.
    3. Swap `array[k]` with `array[l]`.
    4. Reverse the order of the elements from `array[k + 1]` to the end of the array.

    Args:
        array (list): The list of elements to be permuted.

    Returns:
        bool: `True` if the next permutation is generated, `False` if the array is already in the last permutation order.
    """

    # Step 1: Find the largest index k such that array[k] < array[k + 1]
    k = len(array) - 2
    while k >= 0 and array[k] >= array[k + 1]:
        k -= 1

    if k == -1:
        # This is the last permutation
        return False

    # Step 2: Find the largest index l greater than k such that array[k] < array[l]
    l = len(array) - 1
    while array[k] >= array[l]:
        l -= 1

    # Step 3: Swap the value of array[k] with that of array[l]
    array[k], array[l] = array[l], array[k]

    # Step 4: Reverse the sequence from array[k + 1] up to the end
    array[k + 1:] = array[k + 1:][::-1]  # Use slicing to reverse segment

    return True

def gen_next_permute(array):
    """
    Generate all unique permutations of the input array in lexicographical order.

    This function generates and yields all unique permutations of the input array. The permutations
    are generated in lexicographical order, starting with the smallest permutation (in sorted order).

    Args:
        array (list): The list of elements to permute.

    Yields:
        list: The next permutation of the input list.

    Example:
        >>> for perm in gen_next_permute([1, 2, 3]):
        >>>     print(perm)
        [1, 2, 3]
        [1, 3, 2]
        [2, 1, 3]
        [2, 3, 1]
        [3, 1, 2]
        [3, 2, 1]
    """
    #array.sort()            # Start with the smallest lexicographical permutation
    yield array[:]          # Yield the first permutation
    while next_unique_permute(array):
        yield array[:]

def num_bond_permutes(base):
    n = len(base)
    counts = Counter(base)
    numerator = math.factorial(n)
    denominator = 1
    for count in counts.values():
        denominator *= math.factorial(count)
    num_permutations = numerator // denominator
    return num_permutations

# process a single permute #
def process_permute(bond_permute):
    """
    bond_site_dict_2 = {0: (1,0), 1: (1,9), 2: (3,2), 3: (3,9),
                    4: (5,4), 5: (5,8), 6: (7,6), 7: (7,8)}
    """
    Z_sites = np.array([bond_site_dict[bond][0] for bond in bond_permute])
    A_sites = np.array([bond_site_dict[bond][1] for bond in bond_permute])
    A_sites_combi = gen_A_sites_combi(A_sites)

    unique_count = len(np.unique(A_sites))
    if Trie_true:
        Z_fixed_combi_sum_Trie = compute_combi_coeff_Trie(bond_permute, A_sites_combi, Z_sites, A_sites)
        permute_sum = Z_fixed_combi_sum_Trie * 4 * (2 ** unique_count)
    else:
        Z_fixed_combi_sum_classic = compute_combi_coeff_classic(bond_permute, A_sites_combi, Z_sites, A_sites)
        permute_sum = Z_fixed_combi_sum_classic * 4 * (2 ** unique_count)

    return bond_permute, permute_sum

# process a single permute - support #
def gen_A_sites_combi(A_sites):

    # A sites
    A_site_counts = Counter(A_sites)
    A_combis_array = [A_sites]
    for site, count in A_site_counts.items():
        pos_to_replace = np.where(A_sites == site)[0]
        A_combis_array = site_to_ladder(A_combis_array, pos_to_replace, A_combis_dict[count])

    return A_combis_array


def site_to_ladder(combis_array, pos_to_replace, up_down_combis):
    # Ensure that there are positions to replace and permutations to apply

    if len(pos_to_replace) == 0 or len(up_down_combis) == 0:

        return combis_array  # If nothing to replace, return the original array

    replaced_combis = []
    for combi in combis_array:
        for perm in up_down_combis:
            for start_index in range(len(pos_to_replace) - len(perm) + 1):
                # Copy the original array
                new_combi = combi.copy()
                # Place the permutation in the positions of k
                new_combi[pos_to_replace[start_index:start_index + len(perm)]] = perm
                # Append the new array to the list
            replaced_combis.append(new_combi)
    return replaced_combis



def gen_OHS_sign(bond_permute, A_site_combi):
    y_bond_indices = np.where(np.isin(bond_permute, y_bonds))[0]   # no. of y-down (first part)

    # Z-site
    Z_site_combi = np.full(12, -1)
    y_bond_up_down_Z = Z_site_combi[y_bond_indices]
    y_up_count_Z = np.sum(y_bond_up_down_Z == 1)
    y_down_count_Z = np.sum(y_bond_up_down_Z == -1)

    # A-site
    y_bond_up_down_A = A_site_combi[y_bond_indices]
    y_up_count_A = np.sum(y_bond_up_down_A == 1)
    y_down_count_A = np.sum(y_bond_up_down_A == -1)

    # num y_up and y_down
    y_up_count = y_up_count_A + y_up_count_Z
    y_down_count = y_down_count_A + y_down_count_Z
    y_total_count = y_up_count + y_down_count

    sign_1 = 1 if y_total_count % 4 == 0 else -1
    sign_2 = -1 if y_down_count % 2 == 0 else 1

    return sign_1 * sign_2

def compute_combi_coeff_Trie(bond_permute, A_site_combis, Z_sites, A_sites):

    # Initialize the Trie

    #energy_factor, spin_factor = 1,1
    trie = Trie((spin_S,) * 10, 1, 1)

    Z_fixed_combi_sum = 0
    
    """
    for A_site_combi in A_site_combis:
        print(A_site_combi)
    """
    for A_site_combi in A_site_combis:
        tot_energy_factor, tot_spin_factor = trie.insert_and_compute(A_site_combi, Z_sites, A_sites, combi_propagate)

        # compute OHS_sign
        OHS_sign = gen_OHS_sign(bond_permute, A_site_combi)

        # compute combi_coeff
        # minus sign for EHS_sign - see v3
        combi_coeff = -(OHS_sign * tot_spin_factor) / tot_energy_factor

        Z_fixed_combi_sum+=combi_coeff

    return Z_fixed_combi_sum

def combi_propagate(spin_config, step, Z_sites, A_sites, A_site_combi):
    spin_config = list(spin_config)

    # extract sites and present spins
    A_site, Z_site = A_sites[step], Z_sites[step]
    A_ini_spin, Z_ini_spin = spin_config[A_site], spin_config[Z_site]
    A_action, Z_action = A_site_combi[step], -1

    # compute spin ladder operator factor (before square root)
    if A_ini_spin == 0.5 and A_action == -1:
        A_spin_factor = 4
    elif A_ini_spin == -0.5 and A_action == 1:
        A_spin_factor = 4
    else:
        A_spin_factor = 3

    if Z_ini_spin == 0.5 and Z_action == -1:
        Z_spin_factor = 4
    elif Z_ini_spin == -0.5 and Z_action == 1:
        Z_spin_factor = 4
    else:
        Z_spin_factor = 3

    # square later
    spin_factor = A_spin_factor * Z_spin_factor

    # update spin sites
    spin_config[Z_site]+=Z_action
    spin_config[A_site]+=A_action

    # compute energy gap (always positive)
    excited_energy = 0.0
    for site in [0,2,4,6,8,9]:
        excited_energy += spin_config[site] * 1.5
    int_z_bond_energy = spin_config[1] * spin_config[7] + spin_config[3] * spin_config[5]
    excited_energy+=int_z_bond_energy
    energy_factor = abs(excited_energy - gs_energy)


    # spin_config is now updated spin_config
    # energy_factor and spin_factor are local, i.e., not cumulative
    return (spin_config, energy_factor, spin_factor)


def compute_combi_coeff_classic(bond_permute, A_sites_combi, Z_sites, A_sites):

    ini_spin_config = (spin_S,) * 10

    for A_site_combi in A_sites_combi:
        cum_energy_factor, cum_spin_factor = 1,1

        spin_config = list(ini_spin_config).copy()

        # propagator
        for step in range(12):
            spin_config, local_energy_factor, local_spin_factor = combi_propagate(spin_config, step, Z_sites, A_sites, A_site_combi)
            if step != 11:
                cum_energy_factor*=local_energy_factor
            cum_spin_factor*=local_spin_factor

        OHS_sign = gen_OHS_sign(bond_permute, A_site_combi)
        combi_coeff = -(OHS_sign * cum_spin_factor) / cum_energy_factor

    return combi_coeff

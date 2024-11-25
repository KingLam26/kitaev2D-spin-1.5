##### import business #####
from utils_analysis import *
import sys
from decimal import Decimal, getcontext
import numpy as np
from math import pi
from scipy.optimize import fsolve

getcontext().prec = 30


##### extract results and preliminary process #####
# process run script input params
runtime_inputs = sys.argv

if len(runtime_inputs) > 1:
    print(f"Received input: spin {sys.argv[1]}")
    spin_S = sys.argv[1]

# read in Koga and Lam results from text file
koga_dict, lam_dict = read_results(f'results_collated.txt')

# extract results for spin_S
koga_coeff, lam_coeff = koga_dict[spin_S], lam_dict[spin_S]


##### analysis 1: compute J_eff_sum and plot #####
# Generate values of t from 0 to pi/4
t_values = np.linspace(0, pi/4, 10000)

# Calculate coupling constants for each t
norm_J_eff_lam = normalize_J_eff(J_eff_sum_lam, lam_coeff, spin_S, t_values)
norm_J_eff_koga = normalize_J_eff(J_eff_sum_koga, koga_coeff, spin_S, t_values)

# Plot the results
plot_koga_lam_single(t_values, norm_J_eff_lam, norm_J_eff_koga, spin_S)
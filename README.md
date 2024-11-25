# Half-Integer-spin-S Kitaev 2D Honeycomb Model

We consider the half-integer-spin-S Kitaev 2D Honeycomb Model defined by 

$$H_0=-J_z \sum_{\langle i, j\rangle_z} S_i^z S_j^z, \quad V=-J_x \sum_{\langle i, j\rangle_x} S_i^x S_j^x-J_y \sum_{\langle i, j\rangle_y} S_i^y S_j^y,$$

in the anisotropic limit. The effective Hamiltonian as per perturbation theory is given by

$$H_{\text {eff }}=\sum_{n=1}^{\infty} P V\left(\frac{\mathbb{I}-P}{E_0-\mathcal{H}_0} V\right)^{n-1} P,$$

where $n$ refers to the order in the expansion of $H_{\text {eff }}$, and $P$ is a projection operator onto the subspace spanned by the ground states of $H_0$. In the case of integer spins, the lowest-order ($4S^{\textrm{th}}$) non-trivial effective Hamitonian is given by 

$$H_{\text {eff }} = J_{\text {eff }} \sum_i\sigma_i^y(\mathrm{L})\sigma_i^z(\mathrm{U})\sigma_i^y(\mathrm{R})\sigma_i^y(\mathrm{D}),$$

where the sum is taken over all $z$ bonds in the 2D honeycomb lattice, and $\sigma_i^z$, $\sigma_i^y$ is a pseudo-half-spin with two degrees-of-freedom. The goal of this repository is to compute the effective coupling constant $J_{\text {eff }}$ for arbitary spin. The set of all allowed perturbative permutation processes contributing non-trivially to $H_{\text {eff }}$ can be organised into a number of diagrams as introduced in [arXiv:1811.05668](
https://doi.org/10.1103/PhysRevB.99.104408). Our code takes a certain diagram and spin-S as input and returns the diagram coefficient given by the sum of all the contributions it contains. As the number of permutation processes to be considered grows by ~$\mathcal{O}(10^5)$ with each spin increment, the problem quickly becomes computationally intractable. Nevertheless, we provide explicit numerical coefficients for spin-3/2, which successfully replicates the result obtained in [arXiv:1811.05668](
https://doi.org/10.1103/PhysRevB.99.104408). With further code optimization (e.g., ``Numba`` and vectorization), we expect spin-5/2 to be tractable with a large cluster.

## File Structure
The project consists of the following files:
1. ``utils.py``: Contains utility functions and helper code for reading in model parameters, generating combinatorial arrays, fine-tuning parallel processing parameters, and saving results.
2. ``core.py``: Implements core functionalities, including generating and processing unique permutations, calculating spin and energy factors, handling overall signs, returns the contributing coefficient to $J_{\mathrm{eff}}$ for each permute.
3. ``main_parallel.py``: Orchestrates the parallel processing framework, including producers and consumers for handling permutations across multiple cores.
4. ``inputs.py``: Specifications on the spin and list of diagrams that one can pick to compute.

## How to Use
1. **Check Input Parameters**: Check model parameters in ``inputs.py``.
2. **Run the Program**: Execute ``main_parallel.py`` from the command line, ``python main_parallel.py bp_1 4``, where ``bp_1`` specifies the label for bp in the input file corresponding to the diagram to be considered, and ``4`` specifies the number of CPU cores for parallel processing.
3. Processed results are saved to an output file named ``run-1-results.txt``.

## Requirements
1. ``Python 3.x``
2. ``NumPy``
3. ``Multiprocessing``

## Sample Output
```
bp_label: bp_1
start time: 2024-08-23 14:14:32
end time: 2024-08-23 14:15:35
run time: 1.06 minutes
cores: 24
total_sum: -91893.99422304098
```

The diagram amplitude along with the relative sign is given by ``total_sum``.

# Numerical experiments for *The tan Θ theorem for stable gyroscopic systems*

MATLAB codes reproducing all tables and figures of Section 4 of the paper
(S. Miodragović, *The tan Θ theorem for stable gyroscopic systems*).

The bounds are evaluated exactly as stated in the paper: the tangent bound (22)
with the constant √2 (Remark 2.2, θ ≤ π/4) and the spectral separation in
1/λ, and the sine bound (32) of Theorem 4.1. In all examples the eigenpairs are
computed from the Hermitian transformed linearization (A_L, I) and
(Ã_L, Ĩ) with `eig` (for the perturbed pair with the option `'chol'`),
or with `eigs` in the matrix-free experiment.

## Requirements

MATLAB R2020a or later. No toolboxes are required. Each script is
self-contained (the wiresaw matrices are generated inside the scripts, so the
NLEVP collection is not needed).

## Scripts

| Script | Produces | Output files |
|---|---|---|
| `Ex1_particle_clean.m` | Example 1: Table 1, Table 2, data for Figure 1, accuracy data | `Ex1_table1.csv`, `Ex1_table2.csv`, `Ex1_accuracy.csv`, `Ex1_bound_ratio.mat` |
| `make_fig1.m` | Figure 1 (a), (b) from `Ex1_bound_ratio.mat` | `fig1a_rev.eps/png`, `fig1b_rev.eps/png` |
| `Ex2_small_gyro_clean.m` | Example 2: Table 3, Figure 2, accuracy data | `Ex2_table3.csv`, `Ex2_accuracy.csv`, `Ex2_bound_ratio.mat`, `Ex2_fig1.eps/png` |
| `Ex3_wiresaw_clean.m` | Example 3: Figure 3, Table 4, accuracy data | `Ex3_fig3_ratio.mat`, `Ex3_table4.csv`, `Ex3_accuracy.csv`, `Ex3_fig1.eps/png` |
| `wiresaw_matrixfree_tangent_large_scale_fft.m` | Example 3, matrix-free experiment: Table 5 (n up to 10^6) | `wiresaw_tangent_paper_table.csv`, `wiresaw_tangent_modes_table.csv`, `wiresaw_tangent_diagnostics.csv` |

Run each script from the folder that contains it (e.g. `run Ex1_particle_clean`).
Run `make_fig1` after `Ex1_particle_clean`.

## Reproducibility

Random perturbations are generated as ε·½(R+Rᵀ)∘X with R uniformly
distributed on [0,1]. The random seeds used for the results in the paper are
fixed in the scripts:

| Example | seed |
|---|---|
| Example 1 | `rng(852683514)` |
| Example 2 | `rng(1)` |
| Example 3 | `rng(852684860)` |
| matrix-free experiment | `rng(1)` |

Replace the fixed seed by `rng('shuffle')` to draw new perturbations.
The folder `results/` contains the CSV output of the runs used in the paper.

## Matrix-free experiment

Products with the dense gyroscopic matrix G of the wiresaw problem are computed
by the FFT using G = 2v (T + H) diag(1,…,n) with T Toeplitz and H Hankel, in
O(n log n) operations. The gyroscopic perturbation is G̃ = S G S with
S = I + ε diag(d_j). The Frobenius norm of L_K^{-1} δG L_K^{-*} is bounded
from above by an exact banded sum plus a rigorous tail estimate
(`bandWidth = 201`). On a laptop (i7-1165G7) the run for n = 10^6 takes
about two minutes.

## Licence

MIT (see `LICENSE`).

# Error Analysis of Generalized Langevin Dynamics

This repository contains MATLAB code to reproduce the numerical experiments in

> Quanjun Lang and Jianfeng Lu, *Error Analysis of Generalized Langevin Dynamics with Approximated Memory Kernels*.

The code implements and tests trajectory-wise error bounds for generalized Langevin equations (GLEs) with perturbed or approximated memory kernels. It covers both first-order (velocity-only) and second-order (position–velocity) stochastic Volterra equations with exponential and subexponential kernels, and compares empirical trajectory discrepancies with the theoretical rates derived in the paper.

---

## Repository structure

Top-level files and folders:

- `addpaths.m`  
  Adds all subfolders of the repository to the MATLAB path.

- `system_info.m`  
  Creates a structure with basic system parameters (dimension, damping, time step, final time, number of trajectories, etc.).

- `update_sysInfo.m`  
  Helper to update or refine `sysInfo` (time grids, derived parameters).

- `new_kernel_info_structure.m`  
  Helper for defining kernel families used in the numerical experiments.

- `Schur_norm.m`  
  Computes the Schur-type norm of the memory kernels used in the analysis.

### Kernels

- `Kernels/kernel_info.m`  
  Defines memory kernel types and their parameters (e.g., exponential, power-law, oscillatory, cut-off, translated/dilated variants).

- `Kernels/plot_kernels.m`  
  Visualization of kernels in different scales (linear, semilog, log–log).

### Lyapunov distance

- `Lyapunov/compute_Lyapunov_distance.m`  
  Computes the Lyapunov-type distance between trajectories.

- `Lyapunov/Lyapunov_info.m`  
  Stores parameters for the Lyapunov function and related constants.

### Potentials (for second-order GLEs)

- `Potentials/potential_info.m`  
  Defines confining potentials and their gradients.

- `Potentials/plot_potential_2d.m`  
  Plotting utility for potentials in two dimensions.

### Trajectory generation and plotting

- `Traj/get_noise.m`  
  Generates shared Brownian noise / random driving terms for trajectories.

- `Traj/generate_data_first_order_highdim.m`  
  Simulator for high-dimensional first-order GLEs with given kernel and system parameters.

- `Traj/generate_data_second_order_highdim.m`  
  Simulator for high-dimensional second-order GLEs (position–velocity formulation).

- `Traj/plot_trajectory_decay_fit.m`  
  Fits and plots the decay of Lyapunov distances or mean-square errors.

- `Traj/plot_velocity_highdim.m`, `Traj/plot_trajectory_highdim.m`  
  Visualization of sample trajectories and velocities.

### Section 5.1 – Differential–integral inequality

Folder: `Sec_5_1_diff_int_ineq/`

- `first_order_loop_exp.m`  
  Experiment loop for the exponential-kernel case.

- `first_order_loop_powerlaw.m`  
  Experiment loop for the power-law kernel case.

- `exp_loop.mat`, `powerlaw_loop.mat`, `powerlaw_loop_large.mat`  
  Saved data for the corresponding experiments.

- `power_law_decay_ex.pdf`  
  Example figure for the power-law decay.

### Section 5.2 – First-order GLE with exponential kernels

Folder: `Sec_5_2_first_order_GLE/`

Representative files:

- `main_first_order_exponential.m`  
  Main script to reproduce the first-order GLE experiment with exponential kernels.

- `main_first_order_exp.m`  
  Alternative entry script (same section, slightly different setup).

- `main_first_order_powerlaw_old.m`  
  Older script for power-law experiments (kept for reference).

- `four_est_kernel_loop.m`  
  Kernel perturbation loop used to compare four types of approximations.

- `ensemble_test_S_c2_C2.m`  
  Ensemble experiments for constants and rates in the error bounds.

- `new_kernel_info.m`, `exp_data_Oct17.mat`, `README_first_order_GLE.txt` (if present)  
  Additional configuration and stored data for this section.

### Section 5.3 – Second-order GLE

Folder: `Sec_5_3_second_order_GLE/`

Representative files:

- `main_second_order_exp_loop.m`  
  Main loop for second-order GLE experiments (various kernel perturbations).

- `ensemble_test_second_order.m`  
  Ensemble experiments for second-order systems.

- `second_order_exp_data.mat`, `*_loop.mat`, `*_data*.mat`  
  Data files generated or used by the scripts.


### Section 5.4 – First-order GLE with power-law kernels

Folder: `Sec_5_4_first_order_GLE/`

Representative files:

- `main_first_order_polwer_law_loop.m`  
  (Typo in name is intentional, kept for backwards compatibility.) Main loop for power-law kernel experiments in first-order GLEs.

- `ensemble_test_first_order.m`  
  Ensemble tests comparing empirical and theoretical decay rates.

- `*_loop.mat`, `*_data*.mat`  
  Data files for this section.


### Theoretical rates

Folder: `Theoretical_rates/`

- `fit_traj_rate.m`  
  Fits decay rates from trajectory data (e.g., via log–linear regression in time).

- `get_theoretical_rate_exponential.m`  
  Computes theoretical exponential decay rates predicted by the theory.

- `get_theoretical_threshold_powerlaw.m`  
  Computes theoretical thresholds and power-law exponents from the analytical bounds.

---

## Requirements

- MATLAB (tested on recent versions such as R2023a/R2024a).
- No nonstandard toolboxes are required beyond base MATLAB.

---

## Getting started

### Clone or download the repository

    git clone https://github.com/QuanjunLang/error_analysis_GLE.git
    cd error_analysis_GLE

### Set up the MATLAB path

In MATLAB, change the working directory to the repository root and run:

    addpaths;   % add all subfolders to the MATLAB path

### Run a core example

Examples:

    % Section 5.1: differential–integral inequality, exponential kernel
    cd('Sec_5_1_diff_int_ineq');
    first_order_loop_exp;

    % Section 5.1: differential–integral inequality, power-law kernel
    first_order_loop_powerlaw;

    % Section 5.2: first-order GLE with exponential kernels
    cd('../Sec_5_2_first_order_GLE');
    main_first_order_exponential;

    % Section 5.3: second-order GLE
    cd('../Sec_5_3_second_order_GLE');
    main_second_order_exp_loop;

    % Section 5.4: first-order GLE with power-law kernels
    cd('../Sec_5_4_first_order_GLE');
    main_first_order_polwer_law_loop;

Each script generates trajectories, computes Lyapunov distances or mean-square errors, and produces plots that correspond to figures in the paper. Some scripts load precomputed `.mat` files; these are included in the corresponding section folders.

---

## Modifying experiments

- Adjust global system parameters in `system_info.m` and `update_sysInfo.m`.
- Change kernel types and perturbations in `Kernels/kernel_info.m` (and `new_kernel_info_structure.m` when used).
- For second-order systems, modify potentials in `Potentials/potential_info.m`.
- Re-run the relevant `main_*.m` script to regenerate trajectories and plots.

---

## Reproducing figures and tables

Rough mapping between paper sections and code:

- **Section 5.1 (differential–integral inequality):**  
  `Sec_5_1_diff_int_ineq/first_order_loop_exp.m` and `first_order_loop_powerlaw.m`.

- **Section 5.2 (first-order GLE with exponential kernels):**  
  `Sec_5_2_first_order_GLE/main_first_order_exponential.m` (and `main_first_order_exp.m`) together with `four_est_kernel_loop.m` and `ensemble_test_S_c2_C2.m`.

- **Section 5.3 (second-order GLE):**  
  `Sec_5_3_second_order_GLE/main_second_order_exp_loop.m` and `ensemble_test_second_order.m`.

- **Section 5.4 (first-order GLE with power-law kernels):**  
  `Sec_5_4_first_order_GLE/main_first_order_polwer_law_loop.m` and `ensemble_test_first_order.m`.

The scripts in `Theoretical_rates/` and `Traj/plot_trajectory_decay_fit.m` are used across sections to fit empirical rates and compare them to the theoretical predictions.

---

## Citation

If you use this code in your research or teaching, please cite:

> Quanjun Lang and Jianfeng Lu, *Error Analysis of Generalized Langevin Dynamics with Approximated Memory Kernels*.

---

## License

This project is released under the MIT License. See the `LICENSE` file for the full text.


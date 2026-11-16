# PU-Sinc-RBF for Elliptic PDEs

MATLAB implementation of **partition-of-unity RBF meshfree collocation** for two-dimensional elliptic partial differential equations, with patchwise LOOCV shape-parameter selection. Sinc, multiquadric (MQ), and Gaussian kernels are supported.

## Files

- `PU_Sinc_RBF_EDP.m` — the main script, including its local helper functions; no additional `.m` files from the interpolation project are required.

The code includes Poisson, Helmholtz, nonlinear Burgers-type, and singular linear test problems on square or Cassini-type domains.

## Run

Place `PU_Sinc_RBF_EDP.m` in the MATLAB working directory and execute:

```matlab
PU_Sinc_RBF_EDP
```

At the top of the script, select `opt.problem`, `opt.domain`, `opt.rbf_list`, and other experiment settings. Results are saved to `EDP_Results.csv` (RMSE, maximum absolute error `LINF`, CPU time, shape parameters, and conditioning diagnostics). Exact solutions are used to **evaluate errors**, not to construct the numerical PDE solution.

**Requirements:** MATLAB, Statistics and Machine Learning Toolbox (`haltonset`), and Signal Processing Toolbox (`sinc`). Global Optimization Toolbox is needed only for the optional `loocv_ga` mode.

## Reference

Rippa, S. (1999). *An algorithm for selecting a good value for the parameter c in radial basis function interpolation*. *Advances in Computational Mathematics*, 11, 193–210. **DOI:** [10.1023/A:1018975909870](https://doi.org/10.1023/A:1018975909870). Foundational reference for LOOCV-based shape-parameter selection; the PDE implementation adapts this idea to local collocation systems.

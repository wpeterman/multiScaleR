# Package index

## Prepare inputs

Build the kernel inputs and declare how each raster covariate is
summarized before optimization.

- [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  : Kernel Scale Preparation
- [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md)
  : Prepare Kernel Inputs from Group-Specific Rasters
- [`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  [`kernel_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  [`landscape_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  [`surface_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  : Define multiScaleR covariate transformations
- [`estimate_multiscale_ram()`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md)
  : Estimate Memory Requirements for Multiscale Optimization

## Optimize scales of effect

Estimate the scale of effect for each covariate and assess the fit.

- [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  : Multiscale optimization
- [`profile_sigma()`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md)
  : Profile Model Fit Across Sigma Parameter Space
- [`diagnostics()`](https://wpeterman.github.io/multiScaleR/reference/diagnostics.md)
  : Retrieve diagnostics from multiScaleR objects

## Summarize and compare models

- [`summary(`*`<multiScaleR>`*`)`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md)
  : Summarize multiScaleR objects
- [`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md)
  : Scale Distance
- [`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md)
  : multiScaleR model selection
- [`bic_tab()`](https://wpeterman.github.io/multiScaleR/reference/bic_tab.md)
  : multiScaleR model selection
- [`print(`*`<multiScaleR>`*`)`](https://wpeterman.github.io/multiScaleR/reference/print.multiScaleR.md)
  : Print method for multiScaleR
- [`print(`*`<multiScaleR_data>`*`)`](https://wpeterman.github.io/multiScaleR/reference/print.multiScaleR_data.md)
  : Print method for multiScaleR_data
- [`print(`*`<summary_multiScaleR>`*`)`](https://wpeterman.github.io/multiScaleR/reference/print.summary_multiScaleR.md)
  : Print method for summary_multiScaleR

## Plot

- [`plot(`*`<multiScaleR>`*`)`](https://wpeterman.github.io/multiScaleR/reference/plot.multiScaleR.md)
  : Plot method for multiScaleR objects
- [`plot(`*`<sigma_profile>`*`)`](https://wpeterman.github.io/multiScaleR/reference/plot.sigma_profile.md)
  : Plot Sigma Profile
- [`plot_kernel()`](https://wpeterman.github.io/multiScaleR/reference/plot_kernel.md)
  : Plot kernel densities
- [`plot_marginal_effects()`](https://wpeterman.github.io/multiScaleR/reference/plot_marginal_effects.md)
  : Plot Marginal Effects from a Fitted Model

## Spatial projection

Apply the optimized kernels to rasters for spatial prediction.

- [`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md)
  : Apply kernel smoothing to raster layers

## Simulation

- [`sim_rast()`](https://wpeterman.github.io/multiScaleR/reference/sim_rast.md)
  : Simulate spatially autocorrelated raster surfaces
- [`sim_dat()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat.md)
  : Simulate data for optimizing scales of effect
- [`sim_dat_unmarked()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_unmarked.md)
  : Simulate data for optimizing scales of effect with \`unmarked\`

## Example data

- [`count_data`](https://wpeterman.github.io/multiScaleR/reference/count_data.md)
  : Example data frame
- [`landscape_counts`](https://wpeterman.github.io/multiScaleR/reference/landscape_counts.md)
  : Example data frame
- [`hab`](https://wpeterman.github.io/multiScaleR/reference/hab.md) :
  Example raster
- [`landscape`](https://wpeterman.github.io/multiScaleR/reference/landscape.md)
  : Simulated raster
- [`pts`](https://wpeterman.github.io/multiScaleR/reference/pts.md) :
  Spatial sample points
- [`surv_pts`](https://wpeterman.github.io/multiScaleR/reference/surv_pts.md)
  : Spatial sample points

## Package

- [`multiScaleR`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-package.md)
  [`multiScaleR-package`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-package.md)
  : multiScaleR

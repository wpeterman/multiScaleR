# Scale and Center Raster Layers Using Model Parameters

Standardizes raster covariates in a \`terra::SpatRaster\` using mean and
standard deviation values extracted from a fitted model \`multiScaleR\`
model.

## Usage

``` r
scale_center_raster(r, multiScaleR, clamp = FALSE, pct_mx = 0)
```

## Arguments

- r:

  A \`terra::SpatRaster\` containing covariate layers to be scaled. All
  layer names must match those used in the \`multiScaleR\` model.

- multiScaleR:

  A \`multiScaleR\` object created using \`kernel_prep\` or
  \`multiScale_optim\`.

- clamp:

  Logical. If \`TRUE\`, scaled values are clamped to the covariate range
  in the model data.

- pct_mx:

  Numeric. If \`clamp\` is \`TRUE\`, this value specifies the amount
  (percentage; positive or negative) by which to expand/contract the
  min/max range when clamping. Default = 0.

## Value

A \`terra::SpatRaster\` with each layer scaled and optionally clamped.

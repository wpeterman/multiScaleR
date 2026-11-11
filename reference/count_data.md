# Example data frame

Example count data to be used for optimizing scales of effect

## Usage

``` r
data(count_data)
```

## Format

A data frame with 75 rows and 2 columns. Data were simulated from a
Poisson distribution with an intercept of 0.5, a \`hab\` effect of 0.75,
and scale of effect (sigma) of 75.

- y:

  –\> Simulated counts at spatial locations

- hab:

  –\> Scaled and centered weighted mean values from the \`hab\` raster
  at each of the \`pts\`

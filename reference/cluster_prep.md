# Identify the originating package of an S3 or S4 model object and load it on a PSOCK cluster

Identify the originating package of an S3 or S4 model object and load it
on a PSOCK cluster

## Usage

``` r
cluster_prep(model, cl)
```

## Arguments

- model:

  An R model object (S3 or S4)

- cl:

  A PSOCK cluster (created with parallel::makeCluster)

## Value

Invisibly returns the package name loaded

## Details

For internal use

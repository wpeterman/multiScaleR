# Extract the Namespace from a Model Call

This function attempts to extract the package namespace from the call
used to fit a model object, assuming the function was called using the
\`pkg::fun()\` syntax.

## Usage

``` r
extract_namespace(x)
```

## Arguments

- x:

  A fitted model object.

## Value

A character string with the namespace (package name), or \`NULL\` if the
namespace cannot be determined.

## Details

For internal use

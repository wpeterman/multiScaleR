# multiScaleR
```{r, echo=FALSE}
knitr::include_graphics("man/figures/logo.png")
```
## An R package to identify the scale of effect of spatial environmental variables in regression analyses.

Functions to simulate data with known scales of effect are included with this package.

Windows users need to install RTools first. Rtools provides a compiler and some helpers to compile code for R in Windows. Download Rtools from here: <https://cran.r-project.org/bin/windows/Rtools/> and select the appropriate Rtools version (4.0 with R 4.x.x)

To install, right click on the ".exe" file and select "Run as administrator". 

To install this package and all supporting packages needed use all functions, execute the following commands in R:

```         
# Install 'remotes' package, if needed
if(!("remotes" %in% list.files(.libPaths()))) {
      install.packages("remotes", repo = "http://cran.rstudio.com", dep = TRUE) 
} 

remotes::install_github("wpeterman/multiScaleR", 
                        build_vignettes = TRUE) # Download package

library(multiScaleR) # Loads package and the other dependencies
```

A vignette to accompany the package is now available with versions \>= 0.2.0. This document walks through all available functions as well provides worked analyses of simulated data.

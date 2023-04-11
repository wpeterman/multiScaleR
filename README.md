# multiScaleR

## An R package to identify the scale of effect of spatial environmental variables in regression analyses.

Functions to simulate data with known scales of effect are included with this package, but require the installation of `NLMR`, which is no longer supported on CRAN.

Windows users need to install RTools first. Rtools provides a compiler and some helpers to compile code for R in Windows. Download Rtools from here: <https://cran.r-project.org/bin/windows/Rtools/> and select the appropriate Rtools version (4.0 with R 4.x.x)

To install, right click on the Rtools40.exe and select "Run as administrator". During the installation make sure to select "Add Rtools to PATH". Otherwise, accept all defaults for everything else.

To install this package and all supporting packages needed use all functions, execute the following commands in R:

```         
# Install 'devtools' package, if needed
if(!("remotes" %in% list.files(.libPaths()))) {
      install.packages("remotes", repo = "http://cran.rstudio.com", dep = TRUE) 
} 

remotes::install_github("cran/RandomFieldsUtils")
remotes::install_github("cran/RandomFields")
remotes::install_github("ropensci/NLMR")

remotes::install_github("wpeterman/multiScaleR", 
                        build_vignettes = FALSE) # Download package

library(multiScaleR) # Loads package and the other dependencies
```

Eventually, a vignette will be developed to expand on the use of the package (stay tuned!).

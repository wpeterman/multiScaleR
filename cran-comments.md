## Test environments

-   local Windows 11, R 4.5.0
-   R-hub (intel, ubuntu-release, linux (R-devel), windows (R-devel), m1-san (R-devel))
-   GitHub Actions (macOS, Windows, Ubuntu)

## R CMD check results

0 errors \| 0 warnings \| 0 notes

## Changes from Previous Submission

This is a re-submission. In this version, I have:

### Addressed CRAN feedback:

-   **Fixed the Description:** Removed "This package". Expanded the description and added a DOI link to a reference.
-   **Wrapped long-running examples in `\donttest{}`**
-   **Added missing `\value` tags:**
    -   print.multiScaleR.Rd: \value
    -   print.multiScaleR_data.Rd: \value
    -   print.summary_multiScaleR.Rd: \value
-   **Removed examples from unexported functions:** `extract_namespace`, `print.multiScaleR_data`, `scale_center_raster`
-   **Have added `verbose` arguments:** Added a logical `verbose` argument to functions so users can suppress console outputs and updates.
-   **Updated `stop(cat(...))`:** Removed `cat` within `stop` function.

## Downstream Dependencies

There are currently no downstream dependencies on CRAN.

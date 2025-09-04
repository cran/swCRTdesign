
<!-- README.md is generated from README.Rmd. Please edit that file -->

# swCRTdesign

<!-- badges: start -->

<!-- badges: end -->

swCRTdesign includes functions for the design and analysis of stepped
wedge cluster randomized trials.

## Installation

You can install the development version of swCRTdesign like so:

``` r
# install.packages("devtools")
devtools::install_github("jphughes9/swCRTdesign")
```

## Example

This is a basic example which shows you how find the power of a SW trial
with 24 clusters and 5 periods

``` r
library(swCRTdesign)
swPwr(swDsn(c(6,6,6,6)), distn="binomial",
n=120, mu0=0.05, mu1=0.035, tau=0.01, eta=0, rho=0, gamma=0, alpha=0.05, retDATA=FALSE)
#> Warning in swPwr(swDsn(c(6, 6, 6, 6)), distn = "binomial", n = 120, mu0 = 0.05,
#> : Cross-sectional design assumed
#> treatment.1 
#>   0.7861896
```

\#You’ll still need to render `README.Rmd` regularly, to keep
`README.md` up-to-date. \#`devtools::build_readme()` is handy for this.

\#In that case, don’t forget to commit and push the resulting figure
files, so they \#display on GitHub and CRAN.


[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/permDep)](https://cran.r-project.org/package=permDep) [![packageversion](https://img.shields.io/badge/Package%20version-1.0.2-orange.svg?style=flat-square)](commits/master) [![Travis-CI Build Status](https://travis-ci.org/stc04003/permDep.svg?branch=master)](https://travis-ci.org/stc04003/permDep) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/stc04003/permDep?branch=master&svg=true)](https://ci.appveyor.com/project/stc04003/permDep) [![Last-changedate](https://img.shields.io/badge/last%20change-2019--04--23-yellowgreen.svg)](/commits/master)

**permDep**
-----------

<!-- README.md is generated from README.Rmd. Please edit that file -->
***permDep*** implements permutation approaches to test for quasi-independence in left-truncated right-censored survival data.

#### Installation

You can install and load **permDep** from CRAN using

``` r
install.packages("permDep")
library(permDep)
```

You can install permDep from github with:

``` r
## install.packages("devtools")
devtools::install_github("stc04003/permDep")
```

#### References:

-   Chiou, S.H., Qian, J., and Betensky, R.A. (2017). Permutation Test for General Dependent Truncation.
-   Tsai, W.Y., (1990). Testing the Assumption of Independence of Truncation Time and Failure Time. 169--177.

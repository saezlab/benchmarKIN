# benchmarkeR: Evaluation of kinase activity inference tools

<!-- badges: start -->
<!-- badges: end -->

## Overview
benchmarkeR provides two complementary benchmarking approaches to evaluate 
kinase activity inference tools. This has already been applied to a number
of computational algorithms and prior knowledge resources [here](https://github.com/saezlab/kinase_benchmark).

To run the benchmarking approaches yourself please follow these tutorials:

- [Perturbation based benchmark](https://github.com/saezlab/benchmarkeR/blob/main/vignettes/perturbBench.Rmd)
- [CPTAC based benchmark](https://github.com/saezlab/benchmarkeR/blob/main/vignettes/cptacBench.Rmd)

## Installation
To install `benchmarkeR` please run:
```
# Install specific versions of dependencies (currently has to be done separate)
remotes::install_version("Matrix", version="1.6-5", repos="http://cran.us.r-project.org")
remotes::install_version("MASS", version="7.3-60", repos="http://cran.us.r-project.org")

# Install the package from GitHub
remotes::install_github("saezlab/benchmarkeR")
```

## Citation
> Preprint out soon

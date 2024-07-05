# benchmarKIN: Evaluation of kinase activity inference tools

<!-- badges: start -->
<!-- badges: end -->

## Overview
benchmarKIN provides two complementary benchmarking approaches to evaluate 
kinase activity inference tools. This has already been applied to a number
of computational algorithms and prior knowledge resources [here](https://github.com/saezlab/kinase_benchmark).

To run the benchmarking approaches yourself please check out the [documentation](https://benchmarkin.readthedocs.io/en/latest/#).

## Installation
To install `benchmarKIN` please run:
```
# Install specific versions of dependencies (currently has to be done separate)
remotes::install_version("Matrix", version="1.6-5", repos="http://cran.us.r-project.org")
remotes::install_version("MASS", version="7.3-60", repos="http://cran.us.r-project.org")

# Install the package from GitHub
remotes::install_github("saezlab/benchmarKIN")
```

## Citation
> Mueller-Dott, Sophia, Eric J. Jaehnig, Khoi Pham Munchic, Wen Jiang, Tomer M. Yaron-Barir, Sara R. Savage, Martin Garrido-Rodriguez, et al. 2024. “Comprehensive Evaluation of Phosphoproteomic-Based Kinase Activity Inference.” bioRxiv. [https://doi.org/10.1101/2024.06.27.601117](https://doi.org/10.1101/2024.06.27.601117).

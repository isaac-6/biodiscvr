
# biodiscvr: Biomarker Discovery Using Composite Value Ratios

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL v3 or
later](https://img.shields.io/badge/License-GPL%20v3%2B-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/isaac-6/biodiscvr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/isaac-6/biodiscvr/actions/workflows/R-CMD-check.yaml)
[![status](https://joss.theoj.org/papers/f143631ed4b6459a399248867b0da166/status.svg)](https://joss.theoj.org/papers/f143631ed4b6459a399248867b0da166)
<!-- badges: end -->

## Overview

`biodiscvr` provides a framework for discovering and evaluating novel or
optimized biomarkers defined as ratios of composite values derived from
feature sets (e.g., regional measurements from imaging data). It is
particularly suited for analyzing longitudinal multi-cohort datasets.

The core functionality utilizes a Genetic Algorithm (GA) to search the
feature space for optimal numerator and denominator combinations based
on biomarker performance metrics calculated using linear mixed-effects
models (Group Separation, Sample Size Estimates).

Notes: - Intended for R version 4.5.1 and above. - Although the
underlying methodology is domain‑agnostic, the current implementation is
tailored to Alzheimer’s disease (AD) research, reflecting the datasets
and use cases that motivated the package. Several preprocessing steps,
expected column names, and default region conventions are AD‑specific.
Users applying biodiscvr to other domains may need to adjust variable
names or modify preprocessing functions. For a guide on how to adapt
this framework to non-AD datasets, see the Contributing (“Adapting to
Other Research Domains”) section below.

## Installation

You can install the development version of `biodiscvr` from
[GitHub](https://github.com/isaac-6/biodiscvr), with or without
vignettes, with:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("isaac-6/biodiscvr", ref = "paper",  dependencies = TRUE, build_vignettes = TRUE)
```

**Note:** You may need Rtools (Windows) or Xcode Command Line Tools
(macOS) installed. The installation will also attempt to install
required packages listed under `Imports` (like `yaml`, `readr`, `dplyr`,
`lme4`, `longpower`, `GA`, `rlang`).

## Getting Started: Workflow Tutorial

The best way to learn how to use the package is through the included
vignette, which walks through the entire workflow step-by-step using
your own data alongside the package’s default configuration files.

After installing the package with build_vignettes = TRUE, you can access
the vignette with:

``` r
# Load the package first if needed
# library(biodiscvr)

# Main workflow vignette
vignette("biodiscvr", package = "biodiscvr")

# Synthetic data example
vignette("biodiscvr_synth", package = "biodiscvr")
```

The vignette covers: \* Setting up your data directories and output
paths. \* Loading data with `load_datasets()`. \* Preprocessing data
with `preprocess_data()`. \* Generating demographic tables with
`create_demographics_table()`. \* Running discovery experiments with
`run_experiments()`. \* Evaluating literature biomarkers with
`evaluate_literature_biomarkers()`.

## Core Functions

The package provides several key functions for the workflow:

- `load_datasets()`: Loads data from structured directories.
- `preprocess_data()`: Validates data, prepares SUV columns based on
  dictionary matching/renaming. Filters data based on minimum entries
  and inclusion criteria.
- `create_demographics_table()`: Generates demographic summary tables.
- `run_experiments()`: Orchestrates single- and multi-cohort discovery
  runs based on an experiments configuration file.
- `evaluate_literature_biomarkers()`: Evaluates predefined biomarkers on
  the processed data.
- `run_ablation()`: Re-evaluates the output of evaluation (or results),
  iteratively removing regions to improve SSE.

Additional functions to perform single experiments \*
`biodiscvr_single()`: Runs the CVR discovery process for a single
dataset (typically called by `run_experiments`). \*
`biodiscvr_multicohort()`: Runs the CVR discovery process optimized
across multiple datasets (typically called by `run_experiments`).

Consult the help page for each function (e.g., `?run_experiments`) for
detailed argument descriptions.

## References

This package builds upon the methodologies described in: -
Llorente-Saguer, I. et al. (2024). *A data-driven framework for
biomarker discovery applied to optimizing modern clinical and
preclinical trials on Alzheimer’s disease*. Brain Communications, Volume
6, Issue 6. [doi link](https://doi.org/10.1093/braincomms/fcae438)

## Citation

If you use `biodiscvr` in your research, please cite it. You can get the
citation information by running: `citation("biodiscvr")` after
installation.

## Contributing

Contributions to `biodiscvr` are welcome! Whether you are reporting a
bug, suggesting a feature, or adapting the framework to a new clinical
domain, your input is valuable.

### Adapting to Other Research Domains

While the default templates are tailored for Alzheimer’s Disease (AD),
the core discovery logic is domain-agnostic. To adapt the package to
another field (e.g., oncology, cardiology, astronomy), follow these
steps:

1.  **Modify Configuration:** Update `inst/files/config.yaml` to reflect
    your domain-specific column names and desired trial parameters
    (power, significance level). After all, the sample size estimate is
    a measure of longitudinal statistical power.
2.  **Custom Preprocessing:** If your data requires unique cleaning
    steps beyond simple renaming, you can modify or extend the logic in
    `R/preprocess.R`.
3.  **Define New Fitness Metrics:** If your clinical goal is not based
    on t-statistics or sample size estimates (e.g., you want to optimize
    for Cox Proportional Hazards or AUC), you can implement a custom
    fitness function in `R/evaluation.R`.
4.  **Discovery Logic:** The heuristic search logic lives in
    `R/biodiscvr.R` and `R/biodiscvr_multicohort` (currently a genetic
    algorithm).

### Development Workflow

1.  **Fork the repository** and create a feature branch.
2.  **Submit a Pull Request** with a clear description of your changes.
3.  Ensure that any new functions are documented using Roxygen2
    comments.
4.  Run `R CMD check` to ensure the package remains stable.

### Reporting Issues

If you find bugs or have questions about implementation, please visit
the [GitHub Issues page](https://github.com/isaac-6/biodiscvr/issues).

### Code of Conduct

Please note that the ‘biodiscvr’ project is released with a Contributor
Code of Conduct (`code_of_conduct.md`). By contributing to this project,
you agree to abide by its terms.

## License

This package is licensed under the GNU General Public License v3.0 or
later (GPL-3.0-or-later). See the [LICENSE](LICENSE.md) file for
details.

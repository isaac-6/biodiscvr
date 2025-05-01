
<!-- README.md is generated from README.Rmd. Please edit that file -->

# biodiscvr: Biomarker Discovery Using Composite Value Ratios

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/biodiscvr)](https://CRAN.R-project.org/package=biodiscvr) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/isaac-6/biodiscvr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/isaac-6/biodiscvr?branch=main) -->
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/isaac-6/biodiscvr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/isaac-6/biodiscvr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

`biodiscvr` provides a framework for discovering and evaluating novel or
optimized biomarkers defined as ratios of composite values derived from
feature sets (e.g., regional measurements from imaging data). It is
particularly suited for analyzing longitudinal multi-cohort datasets.

The core functionality utilizes a Genetic Algorithm (GA) to search the
feature space for optimal numerator and denominator combinations that
optimise specific biomarker performance metrics. These metrics are
calculated using linear mixed-effects models and typically include:

- **Repeatability:** Error percentage (`Rep`), calculated as the
  standard deviation of model residuals (data is log-transformed).
- **Group Separation:** T-statistics comparing biomarker slopes or
  changes between groups (e.g., `SepAB` for amyloid-positivity effect).
- **Statistical Power:** Sample Size Estimates (`SSE`) required to
  detect a certain effect size, calculated using the `longpower`
  package.

The package includes functions to handle the workflow:

1.  Loading data from structured directories (`load_datasets`).
2.  Checking data integrity and preprocessing datasets based on
    inclusion criteria (`preprocess_data`). 3a. Running the CVR
    discovery process for single datasets (`biodiscvr_single`). 3b.
    Running the CVR discovery process optimised for multiple datasets
    (`biodiscvr_multicohort`). 3c. Logging results to CSV files and
    generating convergence plots.

Optional, but useful for reports: 1. Report/log filtering of data (done
in `preprocess_data`) 2. Report demographics (age, sex, education,
amyloid-positivity, unique sites and individuals) for the studies groups
(all, cognitively unimpaired, and cognitively impaired).

## Installation

You can install the development version of `biodiscvr` from
[GitHub](https://github.com/isaac-6/biodiscvr) with:

``` r
# install.packages("devtools") # If you don't have devtools installed
devtools::install_github("isaac-6/biodiscvr")

# Or using the 'remotes' package:
# install.packages("remotes")
# remotes::install_github("isaac-6/biodiscvr")
```

**Note:** You may need Rtools (Windows) or Xcode Command Line Tools
(macOS) installed for package development dependencies. The installation
will also attempt to install packages listed under `Imports` in the
DESCRIPTION file (like `yaml`, `readr`, `dplyr`, `lme4`, `longpower`,
`GA`, `rlang`).

## Basic Usage Example

Here’s a conceptual example demonstrating a typical call to the
single-dataset discovery function. *Note: This example uses `eval=FALSE`
as running it requires specific data setup, configuration files. It
takes about 1 second per GA generation, with a 12th gen Intel CPU,
without running in parallel.*

``` r
# Load the package
# library(biodiscvr)

# --- Prerequisites (Conceptual - These steps require setup) ---

# 0a. A directory with available subfolders per dataset is required.
# 0b. If your regions/features are different from the suggested names (found in "inst/files/data_suv.csv"), you can add an additional column with your names, matching the default ones.
# 0c. Each cohort should have (currently) two files: data.csv (demographics) and data_suv_bi (regional SUV or SUVR values), with matching rows.
# 0d. Inclusion criteria is defined in "inst/files/config.yaml".

# 1. Load datasets (replace path)
# loaded_data_list <- load_datasets(root_path = "path/to/your/data/root")

# 2. Check and prepare data (replace paths, uses config.yaml, dict_suv.csv)
# checked_output <- check_and_prepare_data(
#   loaded_data = loaded_data_list,
#   files_path = "path/to/your/files_folder", # Where dicts and config live
#   config_filename = "config.yaml",
#   log_directory = "path/to/your/log_folder"
# )
# config <- checked_output$config
# data_list <- checked_output$data # Renamed/checked data

# 3. Preprocess data (using config settings)
# preprocessed_data_list <- preprocess_datasets(
#   checked_data = checked_output,
#   verbose = TRUE
# )

# --- Example Call to biodiscvr_single ---
# (Assuming the above steps have produced preprocessed_data_list and config)

# Select data for one dataset
# target_dataset_data <- preprocessed_data_list$ADNI # Example

# Define parameters (often from config or set directly)
# target_group <- "CI" # Group to evaluate fitness for
# output_csv_name <- "discovery_results"
# output_dir <- "results"
# experiment_id <- "run_01_baseline"
# composition_method <- 1 # (Only 1 -arithmetic mean- implemented and tested) 
# If there is demand for it, it can be updated with other compositions, 
# such as a weighted mean, geometric mean, volume-weighted...)

# Define features (example - using all numeric from SUV data)
# id_col <- config$preprocessing$id_column %||% "RID"
# features_to_use <- setdiff(names(target_dataset_data$data_suv_bi), id_col)

# --- Conceptual Call (won't run due to GA and fitness function complexity) ---
```

``` r
# Conceptual call - requires data setup and implemented fitness function
# Assume config and target_dataset_data exist from steps above

discovery_result_df <- biodiscvr_single(
  dataset_data = target_dataset_data,
  dataset_name = "ADNI", # Name for logging
  group = target_group,
  config = config,
  features = features_to_use, # can leave empty if all features from data_suv_bi
  fixed_numerator_regs = NULL, # optional predefinition of fixed regions
  fixed_denominator_regs = NULL, # optional predefinition of fixed regions
  var_composition = composition_method,
  bilateral = TRUE, # Current version only does bilateral region analysis
  experiment_tag = experiment_id,
  output_csv_name = output_csv, # results file name
  output_dir = output_dir, # Set path to save csv results and plot
  save_plot = TRUE, # Control plot saving
  ga_seed = 42 # For reproducibility
  # Add other GA parameters if needed (passed via ...)
)

# Inspect the results data frame (if run succeeded)
# print(discovery_result_df)
```

*Why `eval=FALSE`?* The code above demonstrates the function call
structure but doesn’t execute because running the GA
(`biodiscvr_single`) requires correctly formatted input data,
configuration files, and the specific implementation details within the
internal `.calculate_fitness` function, which are beyond the scope of a
simple README example.

## Features

- **Configurable Workflow:** Control preprocessing, modeling, and GA
  parameters via YAML configuration files.
- **Multi-Cohort Handling:** Load, preprocess, and discover biomarkers
  within multiple datasets stored in a structured directory format.
- **GA-Based Optimization:** Employs a Genetic Algorithm (`GA` package)
  to search for optimal feature combinations (numerators/denominators)
  defining the biomarker ratio, based on a custom fitness function.
- **LME-Based Evaluation:** Uses linear mixed-effects models (`lme4`
  package) internally to calculate key performance metrics
  (Repeatability, Separation, Sample Size Estimates via `longpower`).
- **Group-Specific Analysis:** Optimize and evaluate metrics
  specifically for defined clinical groups (e.g., CU vs. CI).
- **Results Logging:** Append detailed results from discovery runs to a
  structured CSV file.
- **Convergence Plotting:** Optionally save GA convergence plots.

## References

This package builds upon the methodologies described in: -
Llorente-Saguer, I. et al. (2024). *A data-driven framework for
biomarker discovery applied to optimizing modern clinical and
preclinical trials on Alzheimer’s disease*. Brain Communications, Volume
6, Issue 6. [doi link](https://doi.org/10.1093/braincomms/fcae438)

## Citation

If you use `biodiscvr` in your research, please cite it. You can get the
citation information by running:

``` r
citation("biodiscvr")
```

## Contributing

Please note that the ‘biodiscvr’ project is released with a Contributor
Code of Conduct. By contributing to this project, you agree to abide by
its terms.

## License

This package is licensed under the GPL v3 License. See the
[LICENSE](LICENSE.md) file for details.

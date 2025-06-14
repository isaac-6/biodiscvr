% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_ablation.R
\name{run_ablation}
\alias{run_ablation}
\title{Run Iterative Regional Ablation on Discovered Biomarkers}
\usage{
run_ablation(
  discovered_biomarkers_csv_path,
  prepared_data_list,
  config,
  output_ablation_results_csv_path,
  datasets_for_ablation_eval = NULL,
  groups_for_ablation_eval = NULL,
  features_for_ablation = NULL,
  verbose = TRUE
)
}
\arguments{
\item{discovered_biomarkers_csv_path}{Character string. Path to the input CSV
file containing discovered biomarkers. Expected columns include at least:
\code{regs_numerator}, \code{regs_denominator}, \code{var_composition}, and identifiers
like \code{discovery_dataset}, \code{experiment_tag}, \code{group_evaluated}.}

\item{prepared_data_list}{List. The output from \code{preprocess_data()},
containing filtered data for multiple datasets.}

\item{config}{List. The loaded main configuration object from
\code{preprocess_data$config}.}

\item{output_ablation_results_csv_path}{Character string. Full path to save the
ablation results as a \emph{new} CSV file. If the file exists, it will be
overwritten. Directory created if needed.}

\item{datasets_for_ablation_eval}{Character vector. Names of the datasets within
\code{prepared_data_list} to use for evaluating SSE \emph{during the ablation process}.
Defaults to using all datasets present in \code{prepared_data_list}.}

\item{groups_for_ablation_eval}{Character vector. Groups ("CU", "CI") to evaluate
within each dataset \emph{during the ablation process}. Defaults to \code{c("CU", "CI")}.}

\item{features_for_ablation}{Character vector or NULL. A specific list of feature
names to consider for ablation. If NULL (default), common features across
\code{datasets_for_ablation_eval} will be identified and used. The regions from
the input CSV must be a subset of these.}

\item{verbose}{Logical. Print progress messages? Defaults to TRUE.}
}
\value{
A data frame containing the aggregated results of the ablation process
for all processed biomarkers, which is also saved to
\code{output_ablation_results_csv_path}. Each row represents one initial biomarker
and includes its original definition, initial and final aggregated SSE,
the refined regions, and a summary of the ablation process (e.g., number
of regions removed). Returns \code{NULL} if the input CSV cannot be read or no
biomarkers are processed.
}
\description{
Reads a CSV file of discovered biomarkers (typically from \code{run_experiments}),
and for each biomarker, performs an iterative single-region ablation using
\code{region_ablation} to minimize aggregated Sample Size Estimate (SSE).
Results of the ablation process are saved to a new CSV file.
}
\details{
Workflow:
\enumerate{
\item Reads the \code{discovered_biomarkers_csv_path}.
\item Identifies common features across \code{datasets_for_ablation_eval} if
\code{features_for_ablation} is not provided.
\item Iterates through each row (discovered biomarker) in the input CSV.
\item For each biomarker:
a. Validates that its defined numerator and denominator regions are present
within the set of features available for ablation.
b. Calls \code{region_ablation} with the biomarker's definition,
the specified evaluation datasets/groups, and other parameters.
c. Collects key results from the ablation: initial and final aggregated SSE,
initial and final region counts, and number of regions removed.
\item Combines all results into a single data frame.
\item Saves the aggregated ablation results table to \code{output_ablation_results_csv_path}.
}

This function is intended to assess the stability and identify core components
of promising biomarkers found in earlier discovery stages.
}

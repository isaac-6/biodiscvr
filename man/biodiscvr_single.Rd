% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biodiscover.R
\name{biodiscvr_single}
\alias{biodiscvr_single}
\title{Run Single-Dataset Biomarker Discovery using Genetic Algorithm}
\usage{
biodiscvr_single(
  dataset_data,
  dataset_name,
  group,
  config,
  features = NULL,
  fixed_numerator_regs = NULL,
  fixed_denominator_regs = NULL,
  min_bounds = NULL,
  max_bounds = NULL,
  var_composition = 1,
  bilateral = TRUE,
  experiment_tag = NULL,
  output_csv_name = NULL,
  save_plot = FALSE,
  output_dir = NULL,
  ga_seed = NULL,
  ...
)
}
\arguments{
\item{dataset_data}{List containing the prepared 'data' and 'data_suv_bi'
data frames for the target dataset (e.g., \code{prepared_data_list$ADNI}).}

\item{dataset_name}{Character string. The name of the dataset being processed
(e.g., "ADNI"). Used for metadata logging.}

\item{group}{Character string. The group ("CU" or "CI") whose specific data
characteristics will be used for the primary fitness evaluation by the
internal \code{.calculate_fitness} function.}

\item{config}{List. The loaded configuration object, typically from
\code{preprocess_data$config}. Must contain sections like
\code{model_equations}, \code{power_params}, \code{preprocessing}, and \code{genetic_algorithm}.}

\item{features}{Character vector. Names of columns in \code{dataset_data$data_suv_bi}
to be considered as features by the GA. If NULL (default), uses all numeric
columns excluding the ID column specified in \code{config$preprocessing$id_column}.}

\item{fixed_numerator_regs}{fixed numerator regions.}

\item{fixed_denominator_regs}{fixed denominator regions.}

\item{min_bounds}{Numeric vector or single number. Lower bounds for GA parameters.
Matched to the number of variables in the GA encoding (e.g., number of features
if using real-valued weights). Defaults to 0.1 if NULL.}

\item{max_bounds}{Numeric vector or single number. Upper bounds for GA parameters.
Matched to the number of variables in the GA encoding. Defaults to 2.9 if NULL.}

\item{var_composition}{Numeric. A parameter passed to the internal
\code{.calculate_fitness} function, likely indicating how regions/features
are combined (e.g., 0 = volume-weighed, 1 = arithmetic mean, etc.).}

\item{bilateral}{Logical. Metadata flag indicating if the analysis considers
bilateral regions (logged in the output). Defaults to TRUE.}

\item{experiment_tag}{Character string or NULL. An optional tag to identify this
specific run/experiment batch in the output CSV.}

\item{output_csv_name}{Character string or NULL. If provided, the full path
to a CSV file where the single result row for this run will be appended.
Handles header creation automatically. Directory will be created if needed.}

\item{save_plot}{Boolean, to choose if a plot of the search is saved (fitness vs generation)}

\item{output_dir}{Path to the output folder (to save csv, and plots)}

\item{ga_seed}{Numeric or NULL. Seed for the GA's random number generator
for reproducibility. If NULL, a random seed is used.}

\item{...}{Additional arguments passed directly to the underlying GA function
(e.g., \code{GA::ga}).}
}
\value{
A single-row data frame containing the results for the best solution
found by the GA, matching the structure intended for the output CSV.
Columns include metadata, fitness value, final evaluation metrics (Rep, SepAB,
SSE for the evaluated \code{group}), the best biomarker definition (e.g.,
regions), and GA parameters. Returns \code{NULL} if the GA fails, essential
data is missing, or critical steps cannot be completed.
}
\description{
Performs GA optimization to find biomarker definitions (e.g., feature weights
or region selections) for a single dataset, maximizing a custom fitness
function based on internal evaluation metrics.
}
\details{
This function orchestrates the GA optimization for a single dataset:
\enumerate{
\item \strong{Initialization:} Sets up parameters, identifies features, prepares bounds.
\item \strong{Fitness Function:} Defines an internal wrapper around the package's
\code{.calculate_fitness} function. This wrapper handles:
a. Decoding the GA chromosome (solution vector) into a biomarker definition.
\strong{(PLACEHOLDER - requires specific implementation)}
b. Calculating the biomarker \code{value} based on the definition and the
dataset's SUV data. \strong{(PLACEHOLDER - requires specific implementation)}
c. Preparing the necessary data subset for evaluation.
d. Calling \code{.calculate_fitness} with all required arguments.
e. Returning the negative of the fitness value (since \code{.calculate_fitness}
returns "higher is better" and \code{GA::ga} minimizes).
\item \strong{GA Execution:} Runs the genetic algorithm (using \code{GA::ga} by default)
with the wrapper fitness function and parameters derived from the \code{config}
and function arguments.
\item \strong{Result Processing:} Extracts the best solution found by the GA.
\item \strong{Re-evaluation:} Recalculates the biomarker value using the best solution
and calls the internal \code{.feval_group} function to obtain the final set of
evaluation metrics (Rep, SepAB, SSE) for the specified \code{group}.
\item \strong{Output:} Constructs a single-row data frame containing comprehensive
results and metadata.
\item \strong{CSV Logging (Optional):} If \code{output_csv_path} is provided, appends the
result row to the specified CSV file using the internal \code{.append_to_csv} helper.
}

Ensure the internal \code{.calculate_fitness} function exists within the package
and accepts the arguments passed by the wrapper function defined herein.
The PLACEHOLDER sections for decoding the GA solution and calculating the
biomarker value \emph{must} be implemented according to your specific methodology.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twas.R
\name{twas_pipeline}
\alias{twas_pipeline}
\title{Function to perform TWAS analysis for across multiple contexts.
This function peforms TWAS analysis for multiple contexts for imputable genes within an LD region and summarize the twas results.}
\usage{
twas_pipeline(
  twas_weights_data,
  ld_meta_file_path,
  gwas_meta_file,
  region_block,
  rsq_cutoff = 0.01,
  rsq_pval_cutoff = 0.05,
  rsq_option = c("rsq", "adj_rsq"),
  rsq_pval_option = c("pval", "adj_rsq_pval"),
  mr_pval_cutoff = 0.05,
  mr_coverage_column = "cs_coverage_0.95",
  quantile_twas = FALSE,
  output_twas_data = FALSE,
  event_filters = NULL,
  column_file_path = NULL,
  comment_string = "#"
)
}
\arguments{
\item{twas_weights_data}{List of list of twas weights output from generate_twas_db function.}

\item{region_block}{A string with LD region informaiton of chromosome number, star and end position of LD block conneced with "_".}
}
\value{
A list of list containing twas result table and formatted TWAS data compatible with ctwas_sumstats() function.
\itemize{
  \item{twas_table}{ A dataframe of twas results summary is generated for each gene-contexts-method pair of all methods for imputable genes.}
  \item{twas_data}{ A list of list containing formatted TWAS data.}
}
}
\description{
Function to perform TWAS analysis for across multiple contexts.
This function peforms TWAS analysis for multiple contexts for imputable genes within an LD region and summarize the twas results.
}

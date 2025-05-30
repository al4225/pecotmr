% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/encoloc.R
\name{coloc_wrapper}
\alias{coloc_wrapper}
\title{Colocalization Analysis Wrapper}
\usage{
coloc_wrapper(
  xqtl_file,
  gwas_files,
  xqtl_finemapping_obj = NULL,
  xqtl_varname_obj = NULL,
  xqtl_region_obj = NULL,
  gwas_finemapping_obj = NULL,
  gwas_varname_obj = NULL,
  gwas_region_obj = NULL,
  filter_lbf_cs = FALSE,
  filter_lbf_cs_secondary = NULL,
  prior_tol = 1e-09,
  p1 = 1e-04,
  p2 = 1e-04,
  p12 = 5e-06,
  ...
)
}
\arguments{
\item{xqtl_file}{Path to the xQTL RDS file.}

\item{gwas_files}{Vector of paths to GWAS RDS files.}

\item{xqtl_finemapping_obj}{Optional table name in xQTL RDS files (default 'susie_fit').}

\item{xqtl_varname_obj}{Optional table name in xQTL RDS files (default 'susie_fit').}

\item{xqtl_region_obj}{Optional table name in xQTL RDS files (default 'susie_fit').}

\item{gwas_finemapping_obj}{Optional table name in GWAS RDS files (default 'susie_fit').}

\item{gwas_varname_obj}{Optional table name in GWAS RDS files (default 'susie_fit').}

\item{gwas_region_obj}{Optional table name in GWAS RDS files (default 'susie_fit').}

\item{prior_tol}{When the prior variance is estimated, compare the estimated value to \code{prior_tol} at the end of the computation,
and exclude a single effect from PIP computation if the estimated prior variance is smaller than this tolerance value.}

\item{p1, }{p2, and p12 are results from xqtl_enrichment_wrapper (default 'p1=1e-4, p2=1e-4, p12=5e-6', same as coloc.bf_bf).}

\item{region_obj}{Optional table name of region info in susie_twas output filess (default 'region_info').}
}
\value{
A list containing the coloc results and the summarized sets.
}
\description{
This function processes xQTL and multiple GWAS finemapped data files for colocalization analysis.
}
\examples{
xqtl_file <- "xqtl_file.rds"
gwas_files <- c("gwas_file1.rds", "gwas_file2.rds")
result <- coloc_wrapper(xqtl_file, gwas_files, LD_meta_file_path)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmr_functions.R
\name{findTMRs}
\alias{findTMRs}
\title{Find TSS-Proximal Methylation-Controlled Regulatory Sites (TMRs)}
\usage{
findTMRs(
  correlation_list,
  offset_length = 10,
  p_adjust_method = "fdr",
  p_value_threshold = 0.05,
  smoothing_factor = 0.75,
  min_gapwidth = 150,
  min_meth_sites = 5,
  BPPARAM = BiocParallel::bpparam()
)
}
\arguments{
\item{correlation_list}{A list of data.frames with correlation values between methylation sites
and a transcript as returned by calculateMethSiteTranscriptCors.}

\item{offset_length}{Number of methylation sites added upstream and downstream of a central methylation site to form a window, resulting in a window size of 2*offset_length + 1.
Default value is 10.}

\item{p_adjust_method}{The method to use for p-value adjustment. Should be one of the methods in p.adjust.methods. Default is "fdr".}

\item{p_value_threshold}{The p_value cutoff to use (after correcting p-values with p_adjust_method). Default value is 0.05.}

\item{smoothing_factor}{Smoothing factor for exponential moving average. Should be a value between 0 and 1 with higher
values resulting in a greater degree of smoothing. Default is 0.75.}

\item{min_gapwidth}{Merge TMRs with the same direction separated by less than this number of base pairs. Default value is 150.}

\item{min_meth_sites}{Minimum number of methylation sites that TMRs can contain. Default value is 5.}

\item{BPPARAM}{A BiocParallelParam object for parallel processing. Defaults to \code{BiocParallel::bpparam()}.}
}
\value{
A GRanges object with the location of TMRs.
}
\description{
Find TSS-Proximal Methylation-Controlled Regulatory Sites (TMRs)
}

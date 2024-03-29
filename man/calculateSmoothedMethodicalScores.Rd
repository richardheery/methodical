% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmr_functions.R
\name{calculateSmoothedMethodicalScores}
\alias{calculateSmoothedMethodicalScores}
\title{Calculate methodical score and smooth it using a exponential weighted moving average}
\usage{
calculateSmoothedMethodicalScores(
  correlation_df,
  offset_length = 10,
  smoothing_factor = 0.75
)
}
\arguments{
\item{correlation_df}{A data.frame with correlation values between methylation sites and a transcript as returned by calculateMethSiteTranscriptCors.}

\item{offset_length}{Number of methylation sites added upstream and downstream of a central methylation site to form a window, resulting in a window size of 2*offset_length + 1.
Default value is 10.}

\item{smoothing_factor}{Smoothing factor for exponential moving average. Should be a value between 0 and 1 with higher
values resulting in a greater degree of smoothing. Default is 0.75.}
}
\value{
A GRanges object
}
\description{
Calculate methodical score and smooth it using a exponential weighted moving average
}
\examples{

# Load data.frame with CpG methylation-transcription correlation results for TUBB6
data("tubb6_cpg_meth_transcript_cors", package = "methodical")

# Calculate smoothed Methodical scores from correlation values
smoothed_methodical_scores <- methodical::calculateSmoothedMethodicalScores(tubb6_cpg_meth_transcript_cors)

}

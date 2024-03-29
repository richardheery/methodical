% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_meth_rse_accessory_functions.R
\name{.split_meth_array_file}
\alias{.split_meth_array_file}
\title{Split data from a single methylation array files into chunks}
\usage{
.split_meth_array_file(file, column, file_count, parameters)
}
\arguments{
\item{file}{Path to a methylation array file.}

\item{column}{The current grid column being processed.}

\item{file_count}{The number of the file being processed}

\item{parameters}{A list of parameters for processing the bedgraph.}
}
\value{
Invisibly returns NULL.
}
\description{
Split data from a single methylation array files into chunks
}

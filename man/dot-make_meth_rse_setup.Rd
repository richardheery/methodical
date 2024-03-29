% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_meth_rse_accessory_functions.R
\name{.make_meth_rse_setup}
\alias{.make_meth_rse_setup}
\title{Perform setup for makeMethRSEFromBedgraphs or makeMethRSEFromArrayFiles}
\usage{
.make_meth_rse_setup(
  meth_files,
  meth_sites,
  sample_metadata,
  hdf5_dir,
  dataset_name,
  overwrite,
  chunkdim,
  temporary_dir,
  ...
)
}
\arguments{
\item{meth_files}{A vector of paths to files with methylation values.
Automatically detects if meth_files contain a header if every field in the first line is a character.}

\item{meth_sites}{A GRanges object with the locations of the methylation sites of interest. Any regions in meth_files that are not in meth_sites are ignored.}

\item{sample_metadata}{Sample metadata to be used as colData for the RangedSummarizedExperiment.}

\item{hdf5_dir}{Directory to save HDF5 file. Is created if it doesn't exist. HDF5 file is called assays.h5.}

\item{dataset_name}{Name to give data set in HDF5 file.}

\item{overwrite}{TRUE or FALSE indicating whether to allow overwriting if dataset_name already exists in assays.h5.}

\item{chunkdim}{The dimensions of the chunks for the HDF5 file.}

\item{temporary_dir}{Name to give a temporary directory to store intermediate files. A directory with this name cannot already exist.}

\item{...}{Additional arguments to be passed to HDF5Array::HDF5RealizationSink.}
}
\value{
A list describing the setup to be used for makeMethRSEFromBedgraphs or makeMethRSEFromArrayFiles.
}
\description{
Perform setup for makeMethRSEFromBedgraphs or makeMethRSEFromArrayFiles
}

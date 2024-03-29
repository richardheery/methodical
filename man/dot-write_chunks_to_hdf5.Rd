% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_meth_rse_accessory_functions.R
\name{.write_chunks_to_hdf5}
\alias{.write_chunks_to_hdf5}
\title{Write chunks of data to a HDF5 sink}
\usage{
.write_chunks_to_hdf5(temp_chunk_dirs, files_in_chunks, hdf5_sink, hdf5_grid)
}
\arguments{
\item{temp_chunk_dirs}{A vector giving the temporary directory associated with each chunk.}

\item{files_in_chunks}{A list of files associated with each chunk in the order they should be placed.}

\item{hdf5_sink}{A HDF5RealizationSink.}

\item{hdf5_grid}{A RegularArrayGrid.}
}
\value{
Invisibly returns TRUE.
}
\description{
Write chunks of data to a HDF5 sink
}

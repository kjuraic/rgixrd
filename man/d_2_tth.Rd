% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lattice_parameters.R
\name{d_2_tth}
\alias{d_2_tth}
\title{d_2_tth(d, lambda)}
\usage{
d_2_tth(d, lambda = 1.54056)
}
\arguments{
\item{d}{d-spacing in angstroms}

\item{lambda}{X-ray wavelength in angstroms}
}
\value{
tth 2*theta position of XRD peak
}
\description{
calculate position of XRD peak in 2*theta scale for d-spacing and
             wavelength lambda. d-spacing and wavelength shoud be given in same
             units: angstroms or nanometers
}
\examples{
d_2_tth(1.59, 1.54056)
}
\author{
K. Juraić
}

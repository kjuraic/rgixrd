% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lattice_parameters.R
\name{tetragonal_d_2_a}
\alias{tetragonal_d_2_a}
\title{tetragonal_d_2_a(h, k, l, d, c)}
\usage{
tetragonal_d_2_a(h, k, l, d, c)
}
\arguments{
\item{h}{Muller index}

\item{k}{Muller index}

\item{l}{Muller index}

\item{d}{d-spacing}

\item{c}{lattice paramerr}
}
\value{
a lattice parameter of hexagonal crystall structure
}
\description{
calculate lattice parameter a for tetragona crystall structure
             from d-spacing from known Muller indexes (h,k,l) and lattice
             paramter c. For specioal case when Muller index l = 0, a can be
             calculated independently on Muller index c
}
\examples{
tetragonal_d_2_c(1, 0, 0, 1.59, 1)
}
\author{
K. Juraić
}

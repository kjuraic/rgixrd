% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{gixrd_average}
\alias{gixrd_average}
\title{gixrd_average(xrd_dat)}
\usage{
gixrd_average(xrd_dat, scans_to_average = 1:length(xrd_dat))
}
\arguments{
\item{xrd_dat}{list with GIDRD scans data}

\item{scans_to_average}{array of scans ids to average}
}
\value{
data frame with calculated average
}
\description{
calculate average difractogram from list of GIXRD scans read by
             gixrd_read_file function (Maja difractomerer file format).
}
\examples{
\dontrun{gixrd_average(gixrd_data)}
}

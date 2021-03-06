% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{read_bruker_raw4}
\alias{read_bruker_raw4}
\title{read Bruker raw file (version 4)}
\usage{
read_bruker_raw4(file_name)
}
\arguments{
\item{file_name}{full path of file with XRD data}
}
\value{
data.frame(tth, intensity)
}
\description{
read Bruker XRD raw file (version 4.0)
}
\examples{
\dontrun{read_bruker_raw4("XRD_pattern.raw")}

}
\author{
K. Juraić
}

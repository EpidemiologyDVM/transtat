% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transph.R
\name{print.transph_summary}
\alias{print.transph_summary}
\title{Print summary of fitted transph model}
\usage{
\method{print}{transph_summary}(creg_sum, cdigits = 4, pdigits = 3)
}
\arguments{
\item{cdigits}{The minimum number of significant digits to print for 
parameter point and interval estimates.}

\item{pdigits}{The minimum number of significant digits to print for p-
values. This is passed to \code{\link[base]{format.pval}}.}

\item{treg_sum}{An object of class \code{transph_summary}.}
}
\description{
Prints a summary of a fitted \code{transph} model. The number of digits
to print for parameter estimates and p-values can be specified. The 
p-values are formatted using \code{\link[base]{format.pval}}.
}
\author{
Eben Kenah \email{kenah.1@osu.edu}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bucher.R
\name{reformat}
\alias{reformat}
\title{Reformat \code{maicplus_bucher} alike object}
\usage{
reformat(
  x,
  ci_digits = 2,
  pval_digits = 3,
  show_pval = TRUE,
  exponentiate = FALSE
)
}
\arguments{
\item{x}{a list, structured like a \code{maicplus_bucher} object}

\item{ci_digits}{an integer, number of decimal places for point
estimate and derived confidence limits}

\item{pval_digits}{an integer, number of decimal places to display
Z-test p-value}

\item{show_pval}{a logical value, default is TRUE. If FALSE, p-value will not
be output as the second element of the character vector}

\item{exponentiate}{whether the treatment effect and confidence
interval should be exponentiated. This applies to relative
treatment effects. Default is set to false.}
}
\description{
Reformat \code{maicplus_bucher} alike object
}
\keyword{internal}

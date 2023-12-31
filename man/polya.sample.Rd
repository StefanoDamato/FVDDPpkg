% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyasample.R
\name{polya.sample}
\alias{polya.sample}
\title{Sampling via Polya Urn scheme}
\usage{
polya.sample(n, theta, v = c(), sampling.f)
}
\arguments{
\item{n}{The amount of samples to be drawn.}

\item{theta}{The intensity, in the sense of Bayesian Statistics}

\item{v}{A vector of values, considered to be already drawn from the Polya scheme.}

\item{sampling.f}{A function to sample new values. Its unique argument must express the number of values to draw.}
}
\value{
A vector containing n values extracted.
}
\description{
Sampling via Polya Urn scheme
}
\examples{
polya.sample(10, 2, c(0,1), function(x) rbeta(x,1,1))

}

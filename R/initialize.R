#' Initialize Fleming-Viot dependent Dirichlet Processes by setting hyperparameters
#'
#' @param theta The intensity of the centering measure, in the sense of Bayesian
#' Nonparametrics.
#' @param sampling.f A function to sample from the centering. Its unique argument
#' must be the amount of values to be drawn.
#' @param density.f A function to compute the value of the density function or
#' mass function of the centering. It has to be consistent with `sampling.f`.
#' @param atomic A boolean value stating whether the centering is atomic or not.
#'
#' @return A list containing the input (renamed as `theta`, `P0.sample`,
#' `P0.density`, and `is.atomic`) and three empty slots that will contain the
#' information once the FVDDP is updated with data. In particular, they are:
#' * `y.star`: a vector of unique values
#' * `M`: a matrix of multiplicities, represented as row vectors
#' * `w`: a vector of weights associated to each row of the matrix of multiplicities.
#' Such list repesents a n object of the `fvddp` class.
#' @export
#'
#' @examples
#' #initiization with an atomic measure (Pois(3))
#' initialize(theta = 1, sampling.f = function(x) rpois(x, 3),
#'            density.f = function(x) dpois(x, 3), atomic = TRUE)
#'
#' #initialization with a non-atomic measure (N(-1, 3))
#' initialize(theta = 3, sampling.f = function(x) rnorm(x, -1, 3),
#'            density.f = function(x) dnorm(x, -1, 3), atomic = FALSE)
#'
#' @references{
#' \insertRef{PapaspiliopoulosRuggiero2014}{FVDDPpkg}
#'
#' \insertRef{PapaspiliopoulosRuggieroSpanò2016}{FVDDPpkg}
#' }
initialize = function(theta, sampling.f, density.f, atomic){
  fvddp = list('theta' = theta, 'P0.sample' = sampling.f, 'P0.density' = density.f,
              'is.atomic' = atomic, 'y.star' = c(), 'M'=NULL, 'w' = c())

  #assign the correct class and return
  class(fvddp) = 'fvddp'
  return(fvddp)
}

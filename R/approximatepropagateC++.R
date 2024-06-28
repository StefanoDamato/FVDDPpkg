#approximate the propagation
#' Approximate the propagation of a Fleming-Viot latent signal
#'
#' @param fvddp An instance of class  generated via [FVDDPpkg::initialize()].
#' In order to perform the propagation, the FVDDP has to be fed some data using
#' [update()], at least once.
#' @param delta.t The time of the propagation.
#' @param N The amount of samples to be drawn in order to perform the approximation.
#'
#' @return A object of class `fvddp`. Since this function is a Monte-Carlo based
#' approximation of [FVDDPpkg::propagate()], the outputs are similar.
#' @export
#'
#' @examples
#' #a first example
#' FVDDP = initialize(theta = 1, sampling.f = function(x) rpois(x, 3),
#'                    density.f = function(x) dpois(x, 3), atomic = TRUE)
#' FVDDP = update(FVDDP, c(4,5))
#' approx.propagate(FVDDP, 0.2, 10000)
#'
#' #another example; it does not matter wether P0 is atomic or not
#' FVDDP=initialize(theta = 3, function(x) rnorm(x, -1, 3),
#'                  function(x) dnorm(x, -1, 3), atomic = FALSE)
#' FVDDP = update(FVDDP, c(-1.145, 0.553, 0.553, 0.553))
#' approx.propagate(FVDDP, 0.6, 10000)
#'
#' @references{
#'   \insertRef{AscolaniLijoiRuggiero2021}{FVDDPpkg}
#' }
#'
#' @seealso [approx.propagate()] for a (slower) exact computation.
#'
approx.propagate = function(fvddp, delta.t, N) {

  #check the class of the fvddp
  if (!inherits(fvddp, "fvddp")) stop(deparse(substitute(fvddp)), ' not in "fvddp" class')

  #create the vector lambda
  n.max = max(rowSums(fvddp$M))
  lambda = (0:n.max)*(fvddp$theta + 0:n.max -1)/2

  lst = montecarlo_sample_prop_cpp(fvddp$M, delta.t, N, fvddp$w, lambda)

  fvddp$M = lst[[1]]
  colnames(fvddp$M) = fvddp$y.star
  fvddp$w = lst[[2]]/N

  #return the fvddp list
  return(fvddp)
}


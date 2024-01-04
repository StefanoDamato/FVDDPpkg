#' Compare the performance of a Monte-Carlo estimate with respect to the exact result.
#'
#' @param fvddp.exact An instance of class `fvddp` obtained via smoothing
#' ([FVDDPpkg::smooth()]) or propagation ([FVDDPpkg::propagate()]).
#' @param fvddp.approx  An instance of class `fvddp` obtained using the approximating
#' algorithms for smoothing or propagation, with the same input as `fvddp.exact`.
#' @param remove.unmatched Choose whether the weights associated to multiplicities
#' that are in `fvddp.exact` but not in `fvddp.approx` should be removed in the
#' computation (`TRUE`) or considered to be 0 (`FALSE`).
#'
#' @return A vector whose j-th element is the difference (in absolute value) between
#' the weight of the j-th row of the matrix `M` of `fvddp.exact` and the weight
#' of the row of the matrix `M` of `fvddp.approx` equal to it. The length depends
#' on the value of `remove.unmathced`.
#' @export
#'
#' @examples
#' #iniialize the process
#' FVDDP = initialize(3, function(x) rgamma(x, 2,2),
#'                    function(x) dgamma(x, 2,2), FALSE)
#' FVDDP = update(FVDDP, c(rep(abs(rnorm(2,1, 4)), 2), rexp(2, 0.5)))
#' #perform n exact propagation and an approximate one
#' EXACT = propagate(FVDDP, 0.7)
#' APPROX = approx.propagate(FVDDP, 0.7, 10000)
#' #measure the error with this function
#' error.estimate(fvddp.exact = EXACT, fvddp.approx = APPROX, TRUE)
#'
#' #in order to smoot, create and propagate the signal from the past and from the future
#' FVDDP=initialize(3, function(x) rbinom(x, 10, 0.2),
#'                  function(x) dbinom(x, 10, 0.2), TRUE)
#' FVDDP.PAST = update(FVDDP, c(2,3))
#' FVDDP.FUTURE = update(FVDDP, c(4))
#' FVDDP.FUTURE = propagate(FVDDP.FUTURE, 0.5)
#' FVDDP.FUTURE = update(FVDDP.FUTURE, c(1))
#' #compute an exact and an approximate smoothing
#' EXACT = smooth(FVDDP.PAST, FVDDP.FUTURE, 0.4, 0.3, c(1,3))
#' APPROX = approx.smooth(FVDDP.PAST, FVDDP.FUTURE, 0.4, 0.3, c(1,3), 20000)
#' #compute the error again
#' error.estimate(fvddp.exact = EXACT, fvddp.approx = APPROX)
error.estimate = function(fvddp.exact, fvddp.approx, remove.unmatched = F) {

  #check the class of the fvddp
  if (class(fvddp.exact) != 'fvddp') stop(deparse(substitute(fvddp.exact)), ' not in "fvddp" class')
  if (class(fvddp.approx) != 'fvddp') stop(deparse(substitute(fvddp.approx)), ' not in "fvddp" class')

  #check the identity of all rows
  errors = compute_errors_cpp(fvddp.exact$M, fvddp.exact$w,
                              fvddp.approx$M, fvddp.approx$w, remove.unmatched)

  #in this case, remove the indices related to the NAs
  if (remove.unmatched == T) errors = errors[!is.na(errors)]
  return(abs(errors))
}

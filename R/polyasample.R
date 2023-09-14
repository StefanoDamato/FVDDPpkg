
#' Sampling via Polya Urn scheme
#'
#' @param n The amount of samples to be drawn.
#' @param theta The intensity, in the sense of Bayesian Statistics

#' @param v A vector of values, considered to be already drawn from the Polya scheme.
#' @param sampling.f A function to sample new values. Its unique argument must express the number of values to draw.
#'
#' @return A vector containing n values extracted.
#' @export
#'
#' @examples polya.sample(10, 2, c(0,1), function(x) rbeta(x,1,1))
#'
polya.sample = function(n, theta, v=c(), sampling.f){

  N = length(v)
  v = c(v, rep(NA, n))
  for (j in 1:n) {
    v[N+j] = ifelse(runif(1) < theta/(theta + N+j-1), sampling.f(1),
                    sample(v[1:(N+j-1)], 1))
  }
  return(v)
}

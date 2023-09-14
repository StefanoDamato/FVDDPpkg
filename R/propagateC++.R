#' Propagate the Fleming-Viot latent signal in time
#'
#' @param fvddp An instance of class  generated via [FVDDPpkg::initialize()]. In order to
#' perform the propagation, the FVDDP has to be fed some data using
#' [update()], at least once.
#' @param delta.t  The non-negative time of the propagation. If 0, the returned
#' process is the input.
#'
#' @return A list of the same class to the one given as an input (`fvddp`). The
#' amount of rows of the matrix `M`, as well as the vector of weights, `w`, will
#' increase. The hyperparameters will be the same.
#' @export
#'
#' @examples
#' FVDDP = initialize(1, function(x) rpois(x, 3),
#'                    function(x) dpois(x, 3), TRUE)
#' FVDDP = update(FVDDP, c(4,5))
#' propagate(FVDDP, 0.2)
#'
#' FVDDP = initialize(3, function(x) rnorm(x, -1,3),
#'                    function(x) dnorm(x, -1, 3), FALSE)
#' FVDDP = update(FVDDP, c(-1.145, 0.553, 0.553, 0.553))
#' propagate(FVDDP, 0.6)
#'
#' @references{
#' \insertRef{PapaspiliopoulosRuggieroSpanò2016}{FVDDPpkg}
#' }
#'
#' @seealso [approx.propagate()] for a (faster) Monte-Carlo-based analogous.
propagate = function(fvddp, delta.t) {

  #check the class of the fvddp
  if (class(fvddp) != 'fvddp') stop(deparse(substitute(fvddp)), ' not in "fvddp" class')

  #first we exclude two simple cases
  if (delta.t < 0) stop('negative times not accepted')
  if (delta.t==0) return(fvddp)

  #get the old multiplicity matrix
  M = fvddp$M

  #compute the range for all candidate nodes to belong to L(M)
  max.vect = apply(M, 2, max)
  l = list()
  for (j in 1:length(max.vect)) l[[j]] = 0:max.vect[j]
  N = sum(max.vect)

  #generate the matrix LM of such candidates; rename its columns
  LM = as.matrix(expand.grid(l))
  colnames(LM) = fvddp$y.star

  #compute the vector lambda, noting that indexing starts from 0
  lambda = (0:N)*(fvddp$theta + (0:N) -1)/2

  #create a matrix to save the values computed by the function C
  C.matrix = matrix(nrow=N+1, ncol=N)

  #compute the new weights rowwise. NA weights represent a nodes not in L(M)
  w = compute_new_weights_cpp(M, LM, delta.t, fvddp$w, lambda, C.matrix)

  #save the new matrix M, using the not-NA weighted nodes of LM
  fvddp$M = LM[!is.na(w),]
  if (is.null(dim(fvddp$M))) {
    fvddp$M = as.matrix(fvddp$M)
    colnames(fvddp$M) = fvddp$y.star
  }

  #update the weights and return the new fvddp
  fvddp$w = w[!is.na(w)]
  return(fvddp)
}



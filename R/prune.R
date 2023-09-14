#this function allows to remove the components whose weights is under some treshold
#' Reduce the size of Fleming-Viot Dependent Dirichlet Processes
#'
#' @param fvddp An object of class `fvddp`, generated via [initialize()].
#' @param eps The value behold which the weights are removed, with their components
#' of the mixture. `eps` has to be in the interval (0,1).
#'
#' @return An `fvddp` list of smaller size of the input. Precisely, the components
#' whose weight goes below the treshold `eps` will be removed: the vector `w` and
#' the matrix `M` will have a lower amount of rows; if the latter will include all-zero
#' columns, they will be removed.
#' @export
#'
#' @examples
#' #create a large process
#' FVDDP = initialize(3, function(x) rexp(x, 4),
#'                    function(x) dexp(x, 4), FALSE)
#' FVDDP = update(FVDDP, c(rep(rexp(1, 3), 7), rep(rexp(1, 5), 5), rexp(1, 7), 3))
#' FVDDP = propagate(FVDDP, 1)
#' prune(fvddp = FVDDP, eps = 1e-4)
#'
#' @references{
#' \insertRef{AscolaniLijoiRuggiero2023}{FVDDPpkg}
#' }
prune = function(fvddp, eps){

  #check the class of the fvddp
  if (class(fvddp) != 'fvddp') stop(deparse(substitute(fvddp)), ' not in "fvddp" class')

  #check eps
  if (eps > 1) stop('eps must be positive and < 1')

  if (nrow(fvddp$M) > 1) {
    #check all the weights bigger than epsilon
    idx = which(fvddp$w > eps)

    #drop the information associated to the others
    fvddp$w = fvddp$w[idx]/sum(fvddp$w[idx])
    fvddp$M = fvddp$M[idx,]
  }

  #check whether there are unique values that have disappeared
  remove.y = apply(fvddp$M, 2, function(x) all(x==0))

  #if there are all-zero columns
  if (any(remove.y)) {

    #reduce y.star
    fvddp$y.star = fvddp$y.star[!remove.y]

    #handle the case with more the one one multiplicity
    if (nrow(fvddp$M) > 1) {

      #remove them too from M (and y.star)
      fvddp$M = fvddp$M[,!remove.y]
    }

    #in case M has just one row
    else fvddp$M = matrix(fvddp$M[,!remove.y], nrow=1, dimnames = list(NULL,fvddp$y.star))
  }

  #return the pruned fvddp
  return(fvddp)
}

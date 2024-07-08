
#' Use the predictive structure of the FVDDP to sequentially draw values adn update
#'
#' @param fvddp The instance of class `fvddp` the values are drawn from.
#' @param N The amount of values to draw.
#'
#' @return A vector of length `N` of values obtained using the predictive structure.
#' Precisely, after that any observation is drawn (either from the centering measure
#' or from past observations) the input `fvddp` is modified as if the function
#' [FVDDPpkg::update()] is called, with the new value as second argument.
#' @export
#'
#' @examples
#' #create a dumy process and expoit the predictive structure
#' FVDDP = initialize(7, function(x) rbeta(x, 3,3),
#'                    function(x) dgamma(x, 3,3), FALSE)
#' FVDDP = update(FVDDP, rep(0:1, 2))
#' predictive.struct(fvddp = FVDDP, N = 100)
#'
#' @references{
#' \insertRef{AscolaniLijoiRuggiero2021}{FVDDPpkg}
#' }
predictive.struct = function(fvddp, N) {

  #check the class of the fvddp
  if (!inherits(fvddp, "fvddp")) stop(deparse(substitute(fvddp)), ' not in "fvddp" class')

  #initialize an empty vector
  y = rep(NA, N)

  #create some objects to improve notation
  w = fvddp$w
  y.star = fvddp$y.star
  M = fvddp$M

  #compute the row-sums
  row.sums = rowSums(M)

  #initialize a new vector to store new values, ad iterate over it
  for (k in 1:N) {

    #draw a row from the multiplicity matrix and an uniform
    n = as.vector(M[sample(1:nrow(M), 1, prob=w),])
    u = runif(1)

    #with probability theta/(theta + |n|), draw from P0
    if (u < fvddp$theta/(fvddp$theta + sum(n))) y.new = fvddp$P0.sample(1)

    #with probability |n|/(theta + |n|), draw from the empirical measure on n
    else {
      clust = n != 0

      #if there is just one value that can be drawn, i.e. one cluster
      if (sum(clust) == 1) y.new = y.star[clust]

      #otherwise sample
      else y.new = as.numeric(sample(y.star[clust], size=1, prob=n[clust]))
    }


    #find (if it exists) the position of y.new in y.star
    idx = which(y.star == y.new)

    #if it is not present
    if (length(idx) == 0){

      #update the weights dividing by theta + |n|
      w = w/(fvddp$theta + row.sums + k-1)

      #expand y.star and M in order to inclued the new value
      y.star = append(y.star, y.new)
      M = cbind(M, rep(1, nrow(M)))
    }

    #if y.new is in y.star
    else{

      #in case P0 is atomic, the mumerator includes the mass function
      if (fvddp$is.atomic == TRUE) w = w*((fvddp$theta*fvddp$P0.density(y.new) + M[,idx])/
                                        (fvddp$theta + row.sums + k-1))

      #in case P0 is nonatomic, the polya urn scheme may need a theta in the denominator
      else if (fvddp$is.atomic == FALSE) w = w*(ifelse(M[,idx] == 0, fvddp$theta, M[,idx])/
                                              (fvddp$theta + row.sums + k-1))

      #update the multiplicities of the selected type by 1
      M[,idx] = M[,idx]+1
    }

    #add in the place the new value, and normalize w
    y[k] = y.new
    w=w/sum(w)
  }

  #return the vector of new values
  return(y)
}

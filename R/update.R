#update fvddp (when new observations arrive)
#' Update the FVDDP when new observations are collected
#'
#' @param fvddp An object of class `fvddp`; it can be created via [FVDDPpkg::initialize()].
#' @param y.new A vector of new values to update the process.
#'
#' @return An object which is similar to the one given as an input. In particular,
#' the multiplicities of `y.new` will be added to each row of `M`, and the weights
#' `w` will be multiplied times the probability of drawing `y.new` form each row
#' of the matrix `M` according to Polya urn sampling scheme.
#' @export
#'
#' @examples
#' #initialize and propagate a object
#' FVDDP = initialize(1, function(x) rpois(x, 3),
#'                    function(x) dpois(x, 3), TRUE)
#' update(fvddp = FVDDP, y.new = c(4,5))
#'
#' #in this case, update after a propagation to see the diiffent effect of polya urn on the weights
#' FVDDP=initialize(3, function(x) rnorm(x, -1,3),
#'                  function(x) dnorm(x, -1, 3), FALSE)
#' FVDDP = update(FVDDP, c(-1.145, 0.553, 0.553))
#' FVDDP = propagate(FVDDP, 0.6)
#' update(fvddp = FVDDP, y.new = c(0.553, -0.316, -1.145))
#'
#' @references{
#' \insertRef{PapaspiliopoulosRuggieroSpanò2016}{FVDDPpkg}
#' }
update = function(fvddp, y.new) {

  if (!inherits(fvddp, "fvddp")) stop(deparse(substitute(fvddp)), ' not in "fvddp" class')

  #if the model is empty
  if (is.null(fvddp$y.star)) {

    #update all the terms in the fvddp list
    fvddp$y.star = sort(unique(y.new))
    fvddp$M = t(as.matrix(table(y.new), rows=1))
    fvddp$w = c(1)
  }

  #if some observations have already been collected
  else {

    #get a vector of new unique values
    new.y.star = setdiff(unique(y.new), fvddp$y.star)

    #if there are actually new unique values, we should expand M and y.star
    if (length(new.y.star) != 0){

      #append to M some columns of 0s (as many as the new unique values)
      M = (cbind(fvddp$M, matrix(0, nrow=nrow(fvddp$M), ncol=length(new.y.star),
                                 dimnames = list(NULL,new.y.star))))

      #sort the matrix (according to column labels)
      y.star = sort(c(new.y.star, fvddp$y.star))
      M = M[, as.character(y.star)]

      #this if statement just solves a technical dimensionality problem
      if (is.null(dim(M))) {
        M = t(as.matrix(M))
      }

      #update the vector of unique values
      fvddp$y.star = y.star
    }

    #otherwise, there is no need to modify M and y.star
    else {
      M =fvddp$M
    }

    #create a table for the multiplicities of new observations
    y.new.t= table(factor(y.new, levels=fvddp$y.star))

    #compute the new weights, using polya urn for each row of M
    w = apply(M, 1, polya.urn, theta = fvddp$theta, y.new = y.new.t,
              density.f = fvddp$P0.density, atomic = fvddp$is.atomic) * fvddp$w
    fvddp$w = w/sum(w)

    #update the M matrix
    fvddp$M = matrix(apply(M, 1 , function(x) x+ y.new.t), ncol=length(fvddp$y.star), byrow=T)
  }

  #return the process
  return(fvddp)
}


#this function shows how to compute the propabiliy via polya urn sampling scheme
polya.urn = function(y.new, y.tab, theta, density.f, atomic) {

  #in this case ther is noting to compute
  if (all(y.new == 0)) return(1)

  #in this case, we cannot draw the same observation twice from P_0 (a.s.)
  if(atomic == FALSE) {

    #(theta)^(amount of new values)
    p = prod(c(prod(theta*gamma(y.new[(y.new != 0) & (y.tab == 0)])),

               # values already drawn (with their multiplicities)
               prod(gamma((y.new + y.tab)[(y.new != 0) & (y.tab != 0)])
                    /gamma(y.tab[(y.new != 0) & (y.tab != 0)]))))
  }

  #if P_0 is discrete, we can resample from P_0 too
  else if (atomic == TRUE) {

    #use the vectorization property of R's function
    p = prod(gamma(theta*density.f(as.numeric(names(y.new))[y.new != 0])
                   + (y.new + y.tab)[y.new != 0])
             /gamma(theta*density.f(as.numeric(names(y.new))[y.new != 0])
                    + y.tab[y.new != 0]))
  }

  else stop('Type not accepted')

  #the denominator is given by pochammer's symbol
  return(p/prod((theta + sum(y.tab)):(theta + sum(y.tab) + sum(y.new) -1)))
}

#' Compute the smoothing distribution of the Fleming-Viot latent signal
#'
#' @param fvddp.past An instance of class `fvddp`, progressively updated ad propagated
#' with data referring to past times via [FVDDPpkg::update()] and [FVDDPpkg::propagate()]
#' (or its approximate version, [FVDDPpkg::approx.propagate()]).
#' @param fvddp.future Same as `fvddp.past`, but in this case the propagation has
#' been performed with time data from times later than the one to be estimated. Its
#' hyperparameters must be equals to the ones of `fvddp.past`.
#' @param t.past The time between the last collection of data (in the past) and the
#' time at which the smoothing is performed.
#' @param t.future Same as `t.past`, but in this case it is referred to the future.
#' `t.future` is positive too.
#' @param y.new The data collected at the time the smoothing is performed.
#'
#' @return The function returns an instance of class `fvddp` whose hyperparametrs
#' are the same of `fvddp.past` and `fvddp.future`. The values of `y.star`and  `M`
#' are such that to represent the state of the FVDDP signal in the present time,
#' i.e. the one Which is `t.past` away from the time at which `fvddp.past` i
#' estimated, and is `t.future` away from the time at which `fvddp.future` is ,
#' estimated. Since the computation is usually extemely long, one can rely on the
#' Monte-Carlo approximation provided by [FVDDPpkg::approx.smooth()].
#' @export
#'
#' @examples
#' #create wo process and sequentilly update and propagate them
#' FVDDP = initialize(3, function(x) rbinom(x, 10, .2),
#'                    function(x) dbinom(x, 10, .2), TRUE)
#' #for the past
#' FVDDP.PAST = update(FVDDP, c(2,3))
#' #for the future
#' FVDDP.FUTURE = update(FVDDP, c(4))
#' FVDDP.FUTURE = propagate(FVDDP.FUTURE, 0.5)
#' FVDDP.FUTURE = update(FVDDP.FUTURE, c(1))
#' #get a smoothed  FVDDP merging them (with new values too)
#' smooth(fvddp.past = FVDDP.PAST, fvddp.future = FVDDP.FUTURE,
#'        t.past= 0.4, t.future = 0.3, y.new = c(1,3))
#'
#' @references{
#'   \insertRef{AscolaniLijoiRuggiero2023}{FVDDPpkg}
#' }
#'
#' @seealso [approx.smooth()] for a (faster) approximation based on importance
#' sampling.
smooth = function(fvddp.past, fvddp.future, t.past, t.future, y.new){

  #check the class of the arguments
  if (class(fvddp.past) != 'fvddp') stop(deparse(substitute(fvddp.past)),
                                         ' not in "fvddp" class')
  #check the class of the fvddp
  if (class(fvddp.future) != 'fvddp') stop(deparse(substitute(fvddp.future)),
                                           ' not in "fvddp" class')

  #we need the processes (from both past and future) to have the same theta and P0
  if ((fvddp.past$theta != fvddp.future$theta) |
      (fvddp.past$is.atomic != fvddp.past$is.atomic)) stop('Unmatched parameters')

  #initialize a brand new fvddp
  else fvddp = initialize(theta = fvddp.past$theta, sampling.f = fvddp.past$P0.sample,
                          density.f = fvddp.past$P0.density, atomic = fvddp.past$is.atomic)

  #rename the matrices of multiplicities
  M.past = fvddp.past$M
  M.future = fvddp.future$M

  #get the final value of y.star
  y.star = sort(unique(c(fvddp.past$y.star, fvddp.future$y.star, unique(y.new))))

  #check if some values need to be added on y.star on the past or on the future
  new.y.star.past = setdiff(y.star, fvddp.past$y.star)
  new.y.star.future = setdiff(y.star, fvddp.future$y.star)

  #if it is necessary, add some columns of zeros on M.past
  if (length(new.y.star.past) != 0) {
    M.past = cbind(M.past, matrix(0, nrow=nrow(M.past), ncol=length(new.y.star.past),
                                  dimnames = list(NULL,new.y.star.past)))
    M.past = M.past[,order(as.numeric(colnames(M.past)))]

    #modify the type, if it has been coverted to vector
    if (is.null(dim(M.past))) M.past = matrix(M.past, ncol=length(y.star))
  }

  #if it is necessary, add some columns of zeros on M.past
  if (length(new.y.star.future) !=  0){
    M.future = cbind(M.future, matrix(0, nrow=nrow(M.future), ncol=length(new.y.star.future),
                                      dimnames = list(NULL,new.y.star.future)))
    M.future = M.future[,order(as.numeric(colnames(M.future)))]

    #modify the type, if it has been coverted to vector
    if (is.null(dim(M.future))) M.future = matrix(M.future, ncol=length(y.star))
  }

  #write n as a vector of multiplicities
  n  = table(factor(y.new, levels=y.star))

  #save the vectors of weights
  w.past = fvddp.past$w
  w.future = fvddp.future$w

  #save the maximum multiplicites (from past and for future)
  n.past.max = max(rowSums(M.past))
  n.future.max = max(rowSums(M.future))

  #keep these useful information
  N = length(y.star)
  tot = max(n.past.max, n.future.max)
  theta.P0.star = fvddp$theta*fvddp$P0.density(y.star)


  #these objects are necessary to keep in memory some terms
  lambda = (0:tot)*(fvddp$theta + (0:tot) -1)/2
  C.matrix.past = matrix(nrow=n.past.max+1, ncol=n.past.max)
  C.matrix.future = matrix(nrow=n.future.max+1, ncol=n.future.max)
  nonatomic.w.mat = matrix(nrow=n.past.max+1, ncol=n.future.max+1)

  #apply
  lst = compute_M_w_cpp(M.past, M.future, n, t.past, t.future, w.past, w.future, fvddp$is.atomic,
                   lambda, C.matrix.past, C.matrix.future, fvddp$theta, theta.P0.star,
                   nonatomic.w.mat)

  fvddp$y.star = y.star
  fvddp$M = lst$M
  colnames(fvddp$M) = fvddp$y.star
  fvddp$w = lst$w

  #assign the fvddp to the appropriate class
  class(fvddp) = 'fvddp'

  #return the fvddp list
  return(fvddp)
}



#' Show the data contained within the Fleming-Viot Dependent Dirichlet Process
#'
#' @param object An element of class `fvddp`, created via [FVDDPpkg::initialize()].
#' @param ... Optional arguments for `summary` methods.
#' @param rows Specify whether the rows must be printed. Useful in case `M` is large.
#' @param K Specify whether the values of `K`, the amount of clusters for each row,
#' must be printed.
#'
#' @return The function prints a [base::data.frame()] object (that is, of class
#' `"data.frame"`) where every row is a vector of multiplicities (according to
#' the observations as in the names of the columns), with its associated weight.
#' @exportS3Method
#'
#' @examples
#' #iniialize a simple process and show its summary
#' FVDDP = initialize(2, function(x) rgeom(x, .25),
#'                    function(x) dgeom(x, .25), TRUE)
#' FVDDP = update(FVDDP, rpois(4, 2))
#' FVDDP = propagate(FVDDP, 0.5)
#' summary(FVDDP)
summary.fvddp = function(object, ..., rows = FALSE, K = TRUE){

  #check that the argument of the class is correct
  if (!inherits(object, "fvddp")) stop(deparse(substitute(object)), ' not in "fvddp" class')

  #group M and w in a spaced dataframe
  df = as.data.frame(object$M)
  df$'  ' = rep(' ', length(object$w))
  df$weights = object$w

  #in case K has to be printed, add a space to it
  if (K == TRUE) {
    df$' ' = rep(' ', length(object$w))
    df$K = apply(object$M, 1, function(x) sum(x!=0))
  }
  print(df, row.names=rows)
}

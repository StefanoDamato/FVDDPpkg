
#' Print hyperparameters and values from Fleming-Viot Dependent Dirichlet Processes
#'
#' @param x The `fvddp` object to be printed.
#' @param ... Optional arguments for `summary` methods.
#'
#' @return A list of the hyperparameters of the process, i.e. `theta`, `P0.sample`,
#' `Po.density`, and `is.atomic`. Moreover, if the process is still empty, this
#' will be printed; if otherwise it has been updated (via [FVDDPpkg::update()]),
#' then the values in `y.star`, `M` and `w` will be printed.
#' @exportS3Method
#'
#' @examples
#' #a simple example
#' FVDDP = initialize(theta = 1, sampling.f = function(x) rpois(x, 3),
#'                    density.f = function(x) dpois(x, 3), atomic = TRUE)
#' FVDDP = update(FVDDP, c(4,5))
#' print(FVDDP)
#'
#' #in case there are no data
#' FVDDP=initialize(theta = 3, function(x) rnorm(x, -1, 3),
#'                  function(x) dnorm(x, -1, 3), atomic = FALSE)
#' print(FVDDP)
print.fvddp = function(x, ...){

  #print the hyperparameters
  cat('theta:', x$theta, '\n\nP0.sample: ')
  print(x$P0.sample)
  cat('\nP0.density/mass: ')
  print(x$P0.density)
  cat('\nis.atomic:', x$is.atomic, '\n')
  if (is.null(x$y.star)) cat('\nEmpty Process \n')

  #print the values of y.star, M and w
  else{
    colnames(x$M) = x$y.star
    cat('\nUnique values (y.star): \n')
    print(x$y.star)
    cat('\nMultiplicities (M): \n')
    print(x$M)
    cat('\nWeights (w): \n')
    print(x$w)
  }
  return(invisible(x))
}

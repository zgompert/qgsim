#' Generate a bivariate plot
#'
#' This function takes the output of qgsim and some specifiers and generates a bivariate plot
#' @param qgsim_res the results of a qgsim function call.
#' @param populations either a single number indicating a population to plot, or a vector with two values indicating populations to plot against each other.
#' @param vars either a single value or a vector of with two values indicating the z/theta values to plot against each other (can be z1, z2, theta1, or theta2).
#' @param type same as for base R plot() function
#' @param pch same as for base R plot() function
#'
#' @return the output of an r plot() function displaying the required bivariate plot
qgsim_plot_bivar <- function(qgsim_res, populations, vars, type='p', pch=16) {
  # if needed, replicate value in populations
  if (length(populations) == 1) {populations <- rep(populations, 2)}

  # check that requested populations are valid
  npops <- length(qgsim_res$z)
  if (length(populations) != 2) {
    stop("Invalid length of argument 'populations'. Must be 1 or 2.")
  } else if (max(populations) > npops || min(populations) <1) {
    stop(paste("Invalid value in populations. Must be from 1:", npops, sep=""))
  }

  # if needed, replicate value in vars
  if (length(vars) == 1) {vars <- rep(vars, 2)}

  # check that requested vars are valid
  valid_vars <- c('z1', 'z2', 'theta1', 'theta2')
  invalid_vars <- setdiff(vars, valid_vars)
  if (length(invalid_vars) > 0) {
    stop(paste("Error: Invalid var/vars: ", paste(invalid_vars, collapse=", "), ". Must be in: ", paste(valid_vars, collapse = ", "), ".", sep=""))
  } else if (length(vars) != 2) {
    stop(paste("Error: invalid length for vars (", length(vars), "). Must be 1 or 2.", sep=""))
  }

  title_vals<-data.frame(
    z1='Trait value 1',
    z2='Trait value 2',
    theta1='Optimal trait value 1',
    theta2='Optimal trait value 2'
    )

  # generate x-axis title
  x_title<-paste("Pop.", populations[1], "-", title_vals[[vars[1]]])
  # generate y-axis title
  y_title<-paste("Pop.", populations[2], "-", title_vals[[vars[2]]])

  # extract x-values
  x_vals <- switch(
    vars[1],
    "z1"=qgsim_res$z[[populations[1]]][,1],
    "z2"=qgsim_res$z[[populations[1]]][,2],
    "theta1"=qgsim_res$theta[[populations[1]]][,1],
    "theta2"=qgsim_res$theta[[populations[1]]][,2]
  )

  # extract y-values
  y_vals <- switch(
    vars[2],
    "z1"=qgsim_res$z[[populations[2]]][,1],
    "z2"=qgsim_res$z[[populations[2]]][,2],
    "theta1"=qgsim_res$theta[[populations[2]]][,1],
    "theta2"=qgsim_res$theta[[populations[2]]][,2]
  )

  return(plot(x_vals, y_vals, xlab=x_title, ylab=y_title, type=type, pch=pch))
}

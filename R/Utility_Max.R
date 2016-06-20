#' @title Utility Maximization Function
#'
#' @description
#' Utility Maximization.
#'
#' @param DT_bounds hyperparameters lower and upper bounds to limit the search of the acq max
#' @param GP an object of class GP
#' @param acq Acquisition function type to be used
#' @param y_max The current maximum known value of the target utility function
#' @param kappa tunable parameter kappa to balance exploitation against exploration
#' @param eps tunable parameter epsilon to balance exploitation against exploration
#' @return The arg max of the acquisition function
#' @importFrom stats optim
#' @importFrom data.table data.table setnames
#' @importFrom foreach foreach %do%
#' @importFrom magrittr %>%
#' @keywords internal
#' @export

Utility_Max <- function(DT_bounds, GP, acq = "ucb", y_max, kappa, eps) {
  # Try Different Initial Values
  Mat_tries <- Matrix_runif(100,
                            lower = rep(0, length(DT_bounds[, Lower])),
                            upper = rep(1, length(DT_bounds[, Upper])))
  # Negative Utility Function Minimization
  Mat_optim <- foreach(i = 1:nrow(Mat_tries), .combine = "rbind") %do% {
    optim_result <- optim(par = Mat_tries[i,],
                          fn = Utility,
                          GP = GP, acq = acq, y_max = y_max, kappa = kappa, eps = eps,
                          method = "L-BFGS-B",
                          lower = rep(0, length(DT_bounds[, Lower])),
                          upper = rep(1, length(DT_bounds[, Upper])),
                          control = list(maxit = 100,
                                         factr = 5e11))
    c(optim_result$par, optim_result$value)
  } %>%
    data.table(.) %>%
    setnames(., old = names(.), new = c(DT_bounds[, Parameter], "Negetive_Utility"))
  # Return Best set of Parameters
  argmax <- as.numeric(Mat_optim[which.min(Negetive_Utility), DT_bounds[, Parameter], with = FALSE])
  return(argmax)
}
utils::globalVariables(c("Lower", "Upper", "Parameter", "Negetive_Utility"))

#' @title Utility Computing Function
#'
#' @description
#' Computing Utility.
#'
#' @param x_vec a vector of scaled hyperparameters
#' @param GP an object of class GP
#' @param acq Acquisition function type to be used
#' @param y_max The current maximum known value of the target utility function
#' @param kappa tunable parameter kappa to balance exploitation against exploration
#' @param eps tunable parameter epsilon to balance exploitation against exploration
#' @return negative utility to be minimized
#' @importFrom stats pnorm dnorm
#' @importFrom magrittr %>%
#' @keywords internal
#' @export

Utility <- function(x_vec, GP, acq = "ucb", y_max, kappa, eps) {
  # Gaussian Process Prediction
  GP_Pred <- GPfit::predict.GP(object = GP, xnew = matrix(x_vec, nrow = 1))
  GP_Mean <- GP_Pred$Y_hat
  GP_MSE <- GP_Pred$MSE %>% pmax(., 1e-9)
  # Utility Function Type
  if (acq == "ucb") {
    Utility <- GP_Mean + kappa * sqrt(GP_MSE)
  } else if (acq == "ei") {
    z <- (GP_Mean - y_max - eps) / sqrt(GP_MSE)
    Utility <- (GP_Mean - y_max - eps) * pnorm(z) + sqrt(GP_MSE) * dnorm(z)
  } else if (acq == "poi") {
    z <- (GP_Mean - y_max - eps) / sqrt(GP_MSE)
    Utility <- pnorm(z)
  }
  return(-Utility)
}

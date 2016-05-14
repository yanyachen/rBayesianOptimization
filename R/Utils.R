#' @title Matrix runif
#'
#' @description
#' Generates random deviates for multiple hyperparameters in matrix format.
#'
#' @param n number of observations
#' @param lower lower bounds
#' @param upper upper bounds
#' @return a matrix of original hyperparameters
#' @importFrom stats runif
#' @importFrom foreach foreach %do%
#' @importFrom magrittr %>%
#' @keywords internal
#' @export

Matrix_runif <- function(n, lower, upper) {
  foreach(i = seq_along(lower), .combine = "cbind") %do% {
    runif(n, min = lower[i], max = upper[i]) %>%
      pmin(., upper[i] - sqrt(.Machine$double.eps)) %>%
      pmax(., lower[i] + sqrt(.Machine$double.eps))
  } %>%
    matrix(., nrow = n, ncol = length(lower))
}
utils::globalVariables(c("i", "."))


#' @title Matrix MinMax Scaling
#'
#' @description
#' Transforms hyperparameters by scaling each hyperparameter to a given range.
#'
#' @param mat a matrix of original hyperparameters
#' @param lower lower bounds
#' @param upper upper bounds
#' @return a matrix of scaled hyperparameters
#' @importFrom magrittr %>%
#' @keywords internal
#' @export

Min_Max_Scale_Mat <- function(mat, lower, upper) {
  t((t(mat) - lower) / (upper - lower)) %>%
    pmin(., 1 - .Machine$double.eps) %>%
    pmax(., 0 + .Machine$double.eps)
}
utils::globalVariables(".")


#' @title MinMax Inverse Scaling
#'
#' @description
#' Transforms scaled hyperparameters to original range.
#'
#' @param mat a vector of scaled hyperparameters
#' @param lower lower bounds
#' @param upper upper bounds
#' @return a vector of original hyperparameters
#' @importFrom magrittr %>%
#' @keywords internal
#' @export

Min_Max_Inverse_Scale_Vec <- function(vec, lower, upper) {
  (vec * (upper - lower) + lower) %>%
    pmin(., upper - .Machine$double.eps) %>%
    pmax(., lower + .Machine$double.eps)
}

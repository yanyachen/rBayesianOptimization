#' @title Bayesian Optimization
#'
#' @description
#' Bayesian Optimization of Hyperparameters.
#'
#' @param FUN The function to be maximized. This Function should return a named list with 2 components.
#'   The first component "Score" should be the metrics to be maximized, and the second component "Pred" should be
#'   the validation/cross-validation prediction for ensembling/stacking.
#' @param bounds A named list of lower and upper bounds for each hyperparameter.
#'   The names of the list should be identical to the arguments of FUN.
#'   Please use "L" suffix to indicate integer hyperparameter.
#' @param init_points Number of randomly chosen points to sample the
#'   target function before Bayesian Optimization fitting the Gaussian Process.
#' @param n_iter Total number of times the Bayesian Optimization is to repeated.
#' @param acq Acquisition function type to be used. Can be "ucb", "ei" or "poi".
#' \itemize{
#'   \item \code{ucb} GP Upper Confidence Bound
#'   \item \code{ei} Expected Improvement
#'   \item \code{poi} Probability of Improvement
#' }
#' @param kappa tunable parameter kappa to balance exploitation against exploration,
#'   increasing kappa will make the optimized hyperparameters pursuing exploration.
#' @param eps tunable parameter theta to balance exploitation against exploration,
#'   increasing epsilon will make the optimized hyperparameters are more spread out across the whole range.
#' @param verbose Whether or not to print progress.
#' @param ... Other arguments passed on to \code{\link{GP_fit}}.
#' @return a list of Bayesian Optimization result is returned:
#' \itemize{
#'   \item \code{Best_Par} a named vector of the best hyperparameter set found
#'   \item \code{Best_Value} the value of metrics achieved by the best hyperparameter set
#'   \item \code{History} a \code{data.table} of the bayesian optimization history
#'   \item \code{Pred} a \code{data.table} with validation/cross-validation prediction for each round of bayesian optimization history
#' }
#' @references Jasper Snoek, Hugo Larochelle, Ryan P. Adams (2012) \emph{Practical Bayesian Optimization of Machine Learning Algorithms}
#' @examples
#' # Example 1: Optimization
#' ## Set Pred = 0, as placeholder
#' Test_Fun <- function(x) {
#'   list(Score = exp(-(x - 2)^2) + exp(-(x - 6)^2/10) + 1/ (x^2 + 1),
#'        Pred = 0)
#' }
#' ## Set larger init_points and n_iter for better optimization result
#' OPT_Res <- BayesianOptimization(Test_Fun,
#'                                 bounds = list(x = c(1, 3)),
#'                                 init_points = 2, n_iter = 1,
#'                                 acq = "ucb", kappa = 2.576, eps = 0.0,
#'                                 verbose = TRUE)
#' \dontrun{
#' # Example 2: Parameter Tuning
#' library(xgboost)
#' data(agaricus.train, package = "xgboost")
#' dtrain <- xgb.DMatrix(agaricus.train$data,
#'                       label = agaricus.train$label)
#' cv_folds <- KFold(agaricus.train$label, nfolds = 5,
#'                   stratified = TRUE, seed = 0)
#' xgb_cv_bayes <- function(max.depth, min_child_weight, subsample) {
#'   cv <- xgb.cv(params = list(booster = "gbtree", eta = 0.01,
#'                              max_depth = max.depth,
#'                              min_child_weight = min_child_weight,
#'                              subsample = subsample, colsample_bytree = 0.3,
#'                              lambda = 1, alpha = 0,
#'                              objective = "binary:logistic",
#'                              eval_metric = "auc"),
#'                data = dtrain, nround = 100,
#'                folds = cv_folds, prediction = TRUE, showsd = TRUE,
#'                early.stop.round = 5, maximize = TRUE, verbose = 0)
#'   list(Score = cv$dt[, max(test.auc.mean)],
#'        Pred = cv$pred)
#' }
#' OPT_Res <- BayesianOptimization(xgb_cv_bayes,
#'                                 bounds = list(max.depth = c(2L, 6L),
#'                                               min_child_weight = c(1L, 10L),
#'                                               subsample = c(0.5, 0.8)),
#'                                 init_points = 10, n_iter = 20,
#'                                 acq = "ucb", kappa = 2.576, eps = 0.0,
#'                                 verbose = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export

BayesianOptimization <- function(FUN, bounds, init_points, n_iter, acq = "ucb", kappa = 2.576, eps = 0.0, verbose = TRUE, ...) {
  # Preparation
  DT_bounds <- data.table(Parameter = names(bounds),
                          Lower = sapply(bounds, magrittr::extract2, 1),
                          Upper = sapply(bounds, magrittr::extract2, 2),
                          Type = sapply(bounds, class))
  DT_history <- data.table(matrix(-Inf, nrow = init_points + n_iter, ncol = length(bounds) + 2)) %>%
    setnames(., old = names(.), new = c("Round", names(bounds), "Value"))
  Pred_list <- vector(mode = "list", length = init_points + n_iter)
  # Initialization
  for (i in 1:init_points) {
    # Random Sample Point
    Sys.sleep(time = 1)
    set.seed(as.numeric(Sys.time()))
    This_Par <- Matrix_runif(n = 1, lower = DT_bounds[, Lower], upper = DT_bounds[, Upper]) %>%
      as.vector(.) %>%
      magrittr::inset(., DT_bounds[, Type] == "integer", round(magrittr::extract(., DT_bounds[, Type] == "integer"))) %>%
      magrittr::set_names(., DT_bounds[, Parameter])
    # Function Evaluation
    This_Log <- utils::capture.output({
      This_Time <- system.time({
        This_Score_Pred <- do.call(what = FUN, args = as.list(This_Par))
      })
    })
    # Saving History and Prediction
    data.table::set(DT_history, i = as.integer(i),
                    j = names(DT_history),
                    value = as.list(c(Round = i, This_Par, Value = This_Score_Pred$Score)))
    Pred_list[[i]] <- This_Score_Pred$Pred
    # Printing History
    if (verbose == TRUE) {
      paste(c("elapsed", names(DT_history)),
            c(format(This_Time["elapsed"], trim = FALSE, digits = 0, nsmall = 2),
              format(DT_history[i, "Round", with = FALSE], trim = FALSE, digits = 0, nsmall = 0),
              format(DT_history[i, -"Round", with = FALSE], trim = FALSE, digits = 0, nsmall = 4)),
            sep = " = ", collapse = "\t") %>%
        cat(., "\n")
    }
  }
  # Optimization
  for (j in (init_points + 1):(init_points + n_iter)) {
    # Fitting Gaussian Process
    Par_Mat <- Min_Max_Scale_Mat(as.matrix(DT_history[1:(j - 1), DT_bounds[, Parameter], with = FALSE]),
                                 lower = DT_bounds[, Lower],
                                 upper = DT_bounds[, Upper])
    Rounds_Unique <- setdiff(1:(j - 1), which(duplicated(Par_Mat) == TRUE))
    Value_Vec <- DT_history[1:(j - 1), Value]
    GP <- GPfit::GP_fit(X = Par_Mat[Rounds_Unique, ],
                        Y = Value_Vec[Rounds_Unique], ...)
    # Minimizing Negative Utility Function
    Next_Par <- Utility_Max(DT_bounds, GP, acq = acq, y_max = max(DT_history[, Value]), kappa = kappa, eps = eps) %>%
      Min_Max_Inverse_Scale_Vec(., lower = DT_bounds[, Lower], upper = DT_bounds[, Upper]) %>%
      magrittr::inset(., DT_bounds[, Type] == "integer", round(magrittr::extract(., DT_bounds[, Type] == "integer"))) %>%
      magrittr::set_names(., DT_bounds[, Parameter])
    # Function Evaluation
    This_Log <- utils::capture.output({
      This_Time <- system.time({
        Next_Score_Pred <- do.call(what = FUN, args = as.list(Next_Par))
      })
    })
    # Saving History and Prediction
    data.table::set(DT_history, i = as.integer(j),
                    j = names(DT_history),
                    value = as.list(c(Round = j, Next_Par, Value = Next_Score_Pred$Score)))
    Pred_list[[j]] <- Next_Score_Pred$Pred
    # Printing History
    if (verbose == TRUE) {
      paste(c("elapsed", names(DT_history)),
            c(format(This_Time["elapsed"], trim = FALSE, digits = 0, nsmall = 2),
              format(DT_history[j, "Round", with = FALSE], trim = FALSE, digits = 0, nsmall = 0),
              format(DT_history[j, -"Round", with = FALSE], trim = FALSE, digits = 0, nsmall = 4)), sep = " = ", collapse = "\t") %>%
        cat(., "\n")
    }
  }
  # Computing Result
  Best_Par <- as.numeric(DT_history[which.max(Value), DT_bounds[, Parameter], with = FALSE]) %>%
    magrittr::set_names(., DT_bounds[, Parameter])
  Best_Value <- max(DT_history[, Value])
  Pred_DT <- data.table::as.data.table(Pred_list)
  Result <- list(Best_Par = Best_Par,
                 Best_Value = Best_Value,
                 History = DT_history,
                 Pred = Pred_DT)
  # Printing Best
  cat("\n Best Parameters Found: \n")
  paste(names(DT_history),
        c(format(DT_history[which.max(Value), "Round", with = FALSE], trim = FALSE, digits = 0, nsmall = 0),
          format(DT_history[which.max(Value), -"Round", with = FALSE], trim = FALSE, digits = 0, nsmall = 4)),
        sep = " = ", collapse = "\t") %>%
    cat(., "\n")
  # Return
  return(Result)
}
utils::globalVariables(c("Type", "Value"))

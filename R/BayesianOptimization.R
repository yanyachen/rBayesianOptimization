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
#'   All the sample points in init_grid_dt should be in the range of bounds.
#'   Please use "L" suffix to indicate integer hyperparameter.
#' @param init_grid_dt User specified points to sample the target function, should
#'   be a \code{data.frame} or \code{data.table} with identical column names as bounds.
#'   User can add one "Value" column at the end, if target function is pre-sampled.
#' @param init_points Number of randomly chosen points to sample the
#'   target function before Bayesian Optimization fitting the Gaussian Process.
#' @param n_iter Total number of times the Bayesian Optimization is to repeated.
#' @param acq Acquisition function type to be used. Can be "ucb", "ei" or "poi".
#' \itemize{
#'   \item \code{ucb} GP Upper Confidence Bound
#'   \item \code{ei} Expected Improvement
#'   \item \code{poi} Probability of Improvement
#' }
#' @param kappa tunable parameter kappa of GP Upper Confidence Bound, to balance exploitation against exploration,
#'   increasing kappa will make the optimized hyperparameters pursuing exploration.
#' @param eps tunable parameter epsilon of Expected Improvement and Probability of Improvement, to balance exploitation against exploration,
#'   increasing epsilon will make the optimized hyperparameters are more spread out across the whole range.
#' @param kernel Kernel (aka correlation function) for the underlying Gaussian Process. This parameter should be a list
#'   that specifies the type of correlation function along with the smoothness parameter. Popular choices are square exponential (default) or matern 5/2
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
#' xgb_cv_bayes <- function(max_depth, min_child_weight, subsample) {
#'   cv <- xgb.cv(params = list(booster = "gbtree", eta = 0.01,
#'                              max_depth = max_depth,
#'                              min_child_weight = min_child_weight,
#'                              subsample = subsample, colsample_bytree = 0.3,
#'                              lambda = 1, alpha = 0,
#'                              objective = "binary:logistic",
#'                              eval_metric = "auc"),
#'                data = dtrain, nround = 100,
#'                folds = cv_folds, prediction = TRUE, showsd = TRUE,
#'                early_stopping_rounds = 5, maximize = TRUE, verbose = 0)
#'   list(Score = cv$evaluation_log$test_auc_mean[cv$best_iteration],
#'        Pred = cv$pred)
#' }
#' OPT_Res <- BayesianOptimization(xgb_cv_bayes,
#'                                 bounds = list(max_depth = c(2L, 6L),
#'                                               min_child_weight = c(1L, 10L),
#'                                               subsample = c(0.5, 0.8)),
#'                                 init_grid_dt = NULL, init_points = 10, n_iter = 20,
#'                                 acq = "ucb", kappa = 2.576, eps = 0.0,
#'                                 verbose = TRUE)
#' }
#' @importFrom magrittr %>% %T>% extract extract2 inset
#' @importFrom data.table data.table setnames set setDT :=
#' @export

BayesianOptimization <- function(FUN, bounds, init_grid_dt = NULL, init_points = 0, n_iter, acq = "ucb", kappa = 2.576, eps = 0.0, kernel = list(type = "exponential", power = 2), verbose = TRUE, ...) {
  # Preparation
  ## DT_bounds
  DT_bounds <- data.table(Parameter = names(bounds),
                          Lower = sapply(bounds, extract2, 1),
                          Upper = sapply(bounds, extract2, 2),
                          Type = sapply(bounds, class))
  ## init_grid_dt
  setDT(init_grid_dt)
  if (nrow(init_grid_dt) != 0) {
    if (identical(names(init_grid_dt), DT_bounds[, Parameter]) == TRUE) {
      init_grid_dt[, Value := -Inf]
    } else if (identical(names(init_grid_dt), c(DT_bounds[, Parameter], "Value")) == TRUE) {
      paste(nrow(init_grid_dt), "points in hyperparameter space were pre-sampled\n", sep = " ") %>%
        cat(.)
    } else {
      stop("bounds and init_grid_dt should be compatible")
    }
  }
  ## init_points_dt
  init_points_dt <- Matrix_runif(n = init_points, lower = DT_bounds[, Lower], upper = DT_bounds[, Upper]) %>%
    data.table(.) %T>%
    setnames(., old = names(.), new = DT_bounds[, Parameter]) %T>% {
      if (any(DT_bounds[, Type] == "integer")) {
        set(.,
            j = DT_bounds[Type == "integer", Parameter],
            value = round(extract(., j = DT_bounds[Type == "integer", Parameter], with = FALSE)))
      } else {
        .
      }
    } %T>%
    extract(., j = Value:=-Inf)
  ## iter_points_dt
  iter_points_dt <- data.table(matrix(-Inf, nrow = n_iter, ncol = nrow(DT_bounds) + 1)) %>%
    setnames(., old = names(.), new = c(DT_bounds[, Parameter], "Value"))
  ## DT_history
  DT_history <- rbind(init_grid_dt, init_points_dt, iter_points_dt) %>%
    cbind(data.table(Round = 1:nrow(.)), .)
  ## Pred_list
  Pred_list <- vector(mode = "list", length = nrow(DT_history))
  # Initialization
  for (i in 1:(nrow(init_grid_dt) + nrow(init_points_dt))) {
    if (is.infinite(DT_history[i, Value]) == TRUE) {
      This_Par <- DT_history[i, DT_bounds[, Parameter], with = FALSE]
    } else {
      next
    }
    # Function Evaluation
    This_Log <- utils::capture.output({
      This_Time <- system.time({
        This_Score_Pred <- do.call(what = FUN, args = as.list(This_Par))
      })
    })
    # Saving History and Prediction
    data.table::set(DT_history,
                    i = as.integer(i),
                    j = "Value",
                    value = as.list(c(This_Score_Pred$Score)))
    Pred_list[[i]] <- This_Score_Pred$Pred
    # Printing History
    if (verbose == TRUE) {
      paste(c("elapsed", names(DT_history)),
            c(format(This_Time["elapsed"], trim = FALSE, digits = NULL, nsmall = 2),
              format(DT_history[i, "Round", with = FALSE], trim = FALSE, digits = NULL, nsmall = 0),
              format(DT_history[i, -"Round", with = FALSE], trim = FALSE, digits = NULL, nsmall = 4)),
            sep = " = ", collapse = "\t") %>%
        cat(., "\n")
    }
  }
  # Optimization
  for (j in (nrow(init_grid_dt) + nrow(init_points_dt) + 1):nrow(DT_history)) {
    if (nrow(iter_points_dt) == 0) {
      next
    }
    # Fitting Gaussian Process
    Par_Mat <- Min_Max_Scale_Mat(as.matrix(DT_history[1:(j - 1), DT_bounds[, Parameter], with = FALSE]),
                                 lower = DT_bounds[, Lower],
                                 upper = DT_bounds[, Upper])
    Rounds_Unique <- setdiff(1:(j - 1), which(duplicated(Par_Mat) == TRUE))
    Value_Vec <- DT_history[1:(j - 1), Value]
    GP_Log <- utils::capture.output({
      GP <- GPfit::GP_fit(X = Par_Mat[Rounds_Unique, ],
                          Y = Value_Vec[Rounds_Unique],
                          corr = kernel, ...)
    })
    # Minimizing Negative Utility Function
    Next_Par <- Utility_Max(DT_bounds, GP, acq = acq, y_max = max(DT_history[, Value]), kappa = kappa, eps = eps) %>%
      Min_Max_Inverse_Scale_Vec(., lower = DT_bounds[, Lower], upper = DT_bounds[, Upper]) %>%
      magrittr::set_names(., DT_bounds[, Parameter]) %>%
      inset(., DT_bounds[Type == "integer", Parameter], round(extract(., DT_bounds[Type == "integer", Parameter])))
    # Function Evaluation
    Next_Log <- utils::capture.output({
      Next_Time <- system.time({
        Next_Score_Pred <- do.call(what = FUN, args = as.list(Next_Par))
      })
    })
    # Saving History and Prediction
    data.table::set(DT_history,
                    i = as.integer(j),
                    j = c(DT_bounds[, Parameter], "Value"),
                    value = as.list(c(Next_Par, Value = Next_Score_Pred$Score)))
    Pred_list[[j]] <- Next_Score_Pred$Pred
    # Printing History
    if (verbose == TRUE) {
      paste(c("elapsed", names(DT_history)),
            c(format(Next_Time["elapsed"], trim = FALSE, digits = NULL, nsmall = 2),
              format(DT_history[j, "Round", with = FALSE], trim = FALSE, digits = NULL, nsmall = 0),
              format(DT_history[j, -"Round", with = FALSE], trim = FALSE, digits = NULL, nsmall = 4)), sep = " = ", collapse = "\t") %>%
        cat(., "\n")
    }
  }
  # Computing Result
  Best_Par <- as.numeric(DT_history[which.max(Value), DT_bounds[, Parameter], with = FALSE]) %>%
    magrittr::set_names(., DT_bounds[, Parameter])
  Best_Value <- max(DT_history[, Value], na.rm = TRUE)
  Pred_DT <- data.table::as.data.table(Pred_list)
  Result <- list(Best_Par = Best_Par,
                 Best_Value = Best_Value,
                 History = DT_history,
                 Pred = Pred_DT)
  # Printing Best
  cat("\n Best Parameters Found: \n")
  paste(names(DT_history),
        c(format(DT_history[which.max(Value), "Round", with = FALSE], trim = FALSE, digits = NULL, nsmall = 0),
          format(DT_history[which.max(Value), -"Round", with = FALSE], trim = FALSE, digits = NULL, nsmall = 4)),
        sep = " = ", collapse = "\t") %>%
    cat(., "\n")
  # Return
  return(Result)
}
utils::globalVariables(c("Type", "Value"))

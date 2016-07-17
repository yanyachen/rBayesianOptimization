#' @title K-Folds cross validation index generator
#'
#' @description
#' Generates a list of indices for K-Folds Cross-Validation.
#'
#' @param target Samples to split in K folds.
#' @param nfolds Number of folds.
#' @param stratified whether to apply Stratified KFold.
#' @param seed random seed to be used.
#' @return a list of indices for K-Folds Cross-Validation
#' @importFrom stats na.omit
#' @export

KFold <- function(target, nfolds = 10, stratified = FALSE, seed = 0) {
  # Function for generating index for each fold
  nfold_index <- function(index, nfolds, seed) {
    NA_how_many <- ceiling(length(index) / nfolds) * nfolds - length(index)
    set.seed(seed)
    index_mat <- matrix(c(sample(index), rep(NA, NA_how_many)),
                        ncol = nfolds)
    index_list <- list()
    for(i in 1:nfolds) index_list[[i]] <- as.integer(na.omit(index_mat[,i]))
    return(index_list)
  }
  # Function for concatenate corresponding vector in the list
  list_c <- function(list1, list2) {
    stopifnot(length(list1) == length(list2))
    list12 <- list()
    for(i in seq_along(list1)) list12[[i]] <- c(list1[[i]], list2[[i]])
    return(list12)
  }
  # CV-Folds Generating
  if (stratified == FALSE) {
    index_nfold_list <- nfold_index(index = 1:length(target), nfolds = nfolds, seed = seed)
  } else {
    target_table <- table(target)
    target_list <- list()
    for(i in seq_along(names(target_table))) target_list[[i]] <- which(as.character(target) == names(target_table)[i])
    index_list_of_each_target <- lapply(target_list, nfold_index, nfolds = nfolds, seed = seed)
    index_nfold_list <- Reduce(list_c, index_list_of_each_target)
    for (i in seq_along(index_nfold_list)) index_nfold_list[[i]] <- sort(index_nfold_list[[i]])
  }
  return(index_nfold_list)
}

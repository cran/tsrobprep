#' Perform automatic data cleaning of time series data
#'
#' Returns a matrix or a list of matrices with imputed missing values and
#' outliers. The function automatizes the usage of functions
#' \link[tsrobprep]{model_missing_data}, \link[tsrobprep]{detect_outliers} and
#' \link[tsrobprep]{impute_modelled_data}. The function is designed for
#' numerical data only.
#' @param data an input vector, matrix or data frame of dimension nobs x nvars
#' containing missing values; each column is a variable.
#' @param S a number or vector describing the seasonalities (S_1, ..., S_K) in
#' the data, e.g. c(24, 168) if the data consists of 24 observations per day
#' and there is a weekly seasonality in the data.
#' @param tau the quantile(s) of the missing values to be estimated in the
#' quantile regression. Tau accepts all values in (0,1). If NULL, then the
#' weighted lasso regression is performed.
#' @param no.of.last.indices.to.fix a number of observations in the tail of
#' the data to be fixed, by default set to S.
#' @param indices.to.fix indices of the data to be fixed. If NULL, then it is
#' calculated based on the no.of.last.indices.to.fix parameter. Otherwise, the
#' no.of.last.indices.to.fix parameter is ignored.
#' @param model.missing.pars named list containing additional arguments for the
#' \link[tsrobprep]{model_missing_data} function.
#' @param detect.outliers.pars named list containing additional arguments for
#' the \link[tsrobprep]{detect_outliers} function.
#' @details The function calls \link[tsrobprep]{model_missing_data} to clean
#' the data from missing values, \link[tsrobprep]{detect_outliers} to detect
#' outliers, removes them and finally applies again
#' \link[tsrobprep]{model_missing_data} function. For details see the
#' functions' respective help sections.
#' @return A list which contains a matrix or a list of matrices with imputed
#' missing values or outliers, the indices of the data that were modelled, and
#' the given quantile values.
#' @examples
#' \dontrun{
#' autoclean <- auto_data_cleaning(
#'   data = GBload[,-1], S = c(48, 7*48),
#'   no.of.last.indices.to.fix = dim(GBload)[1],
#'   model.missing.pars = list(consider.as.missing = 0, min.val = 0)
#' )
#' autoclean$replaced.indices
#' }
#' @export
#' @seealso \code{\link[tsrobprep]{model_missing_data}},
#' \link[tsrobprep]{detect_outliers}, \link[tsrobprep]{impute_modelled_data}
auto_data_cleaning <- function(
  data, S, tau = NULL, no.of.last.indices.to.fix = S[1],
  indices.to.fix = NULL,
  model.missing.pars = list(), detect.outliers.pars = list()){

  # save the data as matrix
  data.original <- as.matrix(data)
  ###validate the variables' correctness - basic validation
  if(dim(data.original)[1]==0)
    stop("provided data is of length 0") else
      nobs <- dim(data.original)[1]; nvars <- dim(data.original)[2]

  #define the indices to be fixed
  if(is.null(indices.to.fix)) indices.to.fix <-
      seq_len(no.of.last.indices.to.fix) + nobs - no.of.last.indices.to.fix

  # run model_missing_data for the first time - to get rid of missings
  modelled_missings_first <- do.call(tsrobprep::model_missing_data, c(
    list(data = data.original, S = S, tau = 0.5, no.of.last.indices.to.fix = nobs),
    model.missing.pars))

  # impute the modelled values
  data.imputed <- tsrobprep::impute_modelled_data(modelled_missings_first)

  # save the replaced indices
  replaced.indices <- modelled_missings_first$replaced.indices

  # lists for output
  features <- list()
  outlier.probs <- list()

  # for each variable, perform the outlier detection
  for(var.id in seq_len(nvars)){
    # call detect_outliers
    handled_outliers <- do.call(tsrobprep::detect_outliers, c(
      list(data = data.imputed[,var.id], S = S), detect.outliers.pars))

    # if the index is to be fixed, we save it, otherwise it won't be fixed
    # in the final model_missing_data call
    for(i in seq_along(replaced.indices[[var.id]])){
      replaced.indices[[var.id]][[i]] <-
        intersect(replaced.indices[[var.id]][[i]], indices.to.fix)
    }

    # similarly for the outlier case
    replaced.indices[[var.id]][["outliers"]] <- setdiff(
      indices.to.fix[handled_outliers$outlier.pos[indices.to.fix]],
      unlist(replaced.indices[[var.id]]))

    # set the detected outliers to NA
    data.original[replaced.indices[[var.id]][["outliers"]], var.id] <- NA
    # update the feature and outlier probability lists
    features[[var.id]] <- handled_outliers$features
    outlier.probs[[var.id]] <- handled_outliers$outlier.probs
  }

  # final model_missing_data call - after removing the outliers
  modelled_missings <- do.call(tsrobprep::model_missing_data, c(
    list(data = data.original, tau = 0.5, S = S,
         no.of.last.indices.to.fix = no.of.last.indices.to.fix,
         indices.to.fix = indices.to.fix), model.missing.pars))


  # imputing the clean data
  clean.data <- tsrobprep::impute_modelled_data(object = modelled_missings)
  result <- list(clean.data = clean.data,
                 replaced.indices =  replaced.indices,
                 tau = modelled_missings$tau, features = features,
                 outlier.probs = outlier.probs)

  return(result)
}

#' Perform automatic data cleaning of time series data
#'
#' Returns a matrix or a list of matrices with imputed missing values and
#' outliers. The function automatizes the usage of functions
#' \link[tsrobprep]{model_missing_data}, \link[tsrobprep]{handle_outliers} and
#' \link[tsrobprep]{impute_modelled_data}. The function is designed for numerical
#' data only.
#' @param data an input vector, matrix or data frame of dimension nobs x nvars
#' containing missing values; each column is a variable.
#' @param tau the quantile(s) of the missing values to be estimated in the
#' quantile regression. Tau accepts all values in (0,1), the default is 0.5.
#' @param S a number of observations per period, e.g. per day.
#' @param extreg a vector, matrix or data frame of data containing external
#' regressors.
#' @param ... additional arguments for the \link[tsrobprep]{model_missing_data}
#' and \link[tsrobprep]{handle_outliers} functions.
#' @details The function calls \link[tsrobprep]{model_missing_data} and
#' \link[tsrobprep]{handle_outliers} functions simultaneously and combines the
#' results. For details see the functions' respective help sections.
#' @return A matrix or a list of matrices with imputed missing values or
#' outliers.
#' @examples
#' autoclean <- auto_data_cleaning(data = GBload[,-1], S = 48, tau = 0.5,
#'              no.of.last.indices.to.fix = dim(GBload)[1], consider.as.missing = 0,
#'              min.val = 0)
#' @export
#' @seealso \link[tsrobprep]{model_missing_data}, \link[tsrobprep]{handle_outliers},
#' \link[tsrobprep]{impute_modelled_data}
auto_data_cleaning <- function(data, tau = 0.5, S, extreg = NULL, ...){

  data.original <- as.matrix(data)
  ###validate the variables' correctness - basic validation
  if(dim(data.original)[1]==0) stop("provided data is of length 0") else nobs <- dim(data.original)[1]; nvars <- dim(data.original)[2]
  if(!is.numeric(tau) | min(tau) <=0 | max(tau) >= 1) stop("provided tau is incorrect")
  if(!(S == floor(S) & S>0)) stop("provided S is incorrect")
  if(!is.null(extreg) & (is.numeric(extreg) | is.matrix(extreg) | is.data.frame(extreg))) if(nobs != dim(as.matrix(extreg))[1])
    stop("provided data and external regressors are of different number of observations")

  modelled_missings <- tsrobprep::model_missing_data(data = data, tau = tau, S = S, extreg = extreg, ...)

  handled_outliers <- tsrobprep::handle_outliers(data = data, tau = tau, S = S, extreg = extreg, ...)

  full_clean <- modelled_missings

  for(var.id in seq_len(dim(full_clean$original.data)[2])){
    full_clean$replaced.indices[[var.id]] <- c(unlist(modelled_missings$replaced.indices[[var.id]]),
                                               unlist(handled_outliers$replaced.indices[[var.id]]))

    for(tau.no in seq_along(tau)){
      full_clean$replaced.values[[var.id]][unlist(handled_outliers$replaced.indices[[var.id]]), tau.no] <-
        handled_outliers$replaced.values[[var.id]][unlist(handled_outliers$replaced.indices[[var.id]]), tau.no]
    }
  }
  results <- tsrobprep::impute_modelled_data(object = full_clean)
}

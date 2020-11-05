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
#' @param no.of.last.indices.to.fix a number of observations in the tail of
#' the data to be fixed, by default set to S.
#' @param indices.to.fix indices of the data to be fixed. If NULL, then it is
#' calculated based on the no.of.last.indices.to.fix parameter. Otherwise, the
#' no.of.last.indices.to.fix parameter is ignored.
#' @param outlier.lower.cap a number of observations that should fall above
#' (below) the outlier lower cap, used for the calculation of the interquantile
#' range. By default set to 2log(nobs).
#' @param outlier.upper.cap a number of observations that should fall above
#' (below) the outlier upper cap, used for the calculation of the interquantile
#' range. By default set to log(nobs).
#' @param margin a value that indicates the number of interquantile ranges
#' above (below) the quantile upper cap that is still considered not to be an
#' outlier. By default equals outlier.upper.cap.
#' @param max.periods.to.smooth maximum number of S periods to be smoothed in
#' the detection of outliers in the differential series. By default set to 48.
#' @param lags a numeric vector with the lags to use in the autoregression.
#' Negative values are accepted and then also the "future" observations are
#' used for modelling. The default values are constructed under the assumption
#' that S describes daily periodicity.
#' @param extreg a vector, matrix or data frame of data containing external
#' regressors; each column is a variable.
#' @param n.best.extreg a numeric value specifying the maximal number of considered
#' best correlated external regressors (selected in decreasing order). If NULL,
#' then all variables in extreg are used for modelling.
#' @param use.data.as.ext logical specifying whether to use the remaining
#' variables in the data as external regressors or not.
#' @param lag.externals logical specifying whether to lag the external
#' regressors or not. If TRUE, then the algorithm uses the lags specified in
#' parameter lags.
#' @param consider.as.missing a vector of numerical values which are considered
#' as missing in the data.
#' @param whole.period.missing.only if FALSE, then all observations which correspond to
#' the values of consider.as.missing are treated as missings. If TRUE, then only
#' consecutive observations of specified length are considered (length is defined by S).
#' @param min.val a single value or a vector of length nvars providing the
#' minimum possible value of each variable in the data. If a single value, then
#' it applies to all variables. By default set to -Inf.
#' @param max.val a single value or a vector of length nvars providing the
#' maximum possible value of each variable in the data. If a single value, then
#' it applies to all variables. By default set to Inf.
#' @param digits integer indicating the number of decimal places allowed
#' in the data, by default set to 3.
#' @param ... additional arguments for the \link[quantreg]{rq.fit.fnb} algorithm.
#' @details The function calls \link[tsrobprep]{handle_outliers} to detect outliers,
#' removes them and applies \link[tsrobprep]{model_missing_data} function.
#' For details see the functions' respective help sections.
#' @return A matrix or a list of matrices with imputed missing values or
#' outliers.
#' @examples
#' autoclean <- auto_data_cleaning(data = GBload[,-1], S = 48, tau = 0.5,
#'              no.of.last.indices.to.fix = dim(GBload)[1], consider.as.missing = 0,
#'              min.val = 0)
#' @export
#' @seealso \link[tsrobprep]{model_missing_data}, \link[tsrobprep]{handle_outliers},
#' \link[tsrobprep]{impute_modelled_data}
auto_data_cleaning <- function(data, tau = 0.5, S, extreg = NULL, no.of.last.indices.to.fix = S, indices.to.fix = NULL,
                               outlier.lower.cap = 2*log(dim(as.matrix(data))[1]), outlier.upper.cap = log(dim(as.matrix(data))[1]),
                               margin = outlier.upper.cap, max.periods.to.smooth = 48, lags = c(1,2,S,S+1, 7*S, 7*S+1, -1, -2, -S, -S-1, -7*S, -7*S-1),
                               n.best.extreg = NULL, use.data.as.ext = FALSE, lag.externals = FALSE, consider.as.missing = NULL,
                               whole.period.missing.only = FALSE, min.val = -Inf, max.val = Inf, digits = 3, ...){

  data.original <- as.matrix(data)
  ###validate the variables' correctness - basic validation
  if(dim(data.original)[1]==0) stop("provided data is of length 0") else nobs <- dim(data.original)[1]; nvars <- dim(data.original)[2]
  if(!is.numeric(tau) | min(tau) <=0 | max(tau) >= 1) stop("provided tau is incorrect")
  if(!(S == floor(S) & S>0)) stop("provided S is incorrect")
  if(!is.null(extreg) & (is.numeric(extreg) | is.matrix(extreg) | is.data.frame(extreg))) if(nobs != dim(as.matrix(extreg))[1])
    stop("provided data and external regressors are of different number of observations")

  handled_outliers <- tsrobprep::handle_outliers(data = data.original, tau = 0.5, S = S, extreg = extreg, no.of.last.indices.to.fix = no.of.last.indices.to.fix,
                                                 indices.to.fix = indices.to.fix, outlier.lower.cap = outlier.lower.cap,
                                                 outlier.upper.cap = outlier.upper.cap, margin = margin, max.periods.to.smooth = max.periods.to.smooth,
                                                 use.data.as.ext = use.data.as.ext, min.val = min.val, max.val = max.val, ...)
  for(i in seq_len(nvars)){
    data.original[unlist(handled_outliers$replaced.indices[[i]]), i] <- NA
  }

  modelled_missings <- tsrobprep::model_missing_data(data = data.original, tau = tau, S = S, no.of.last.indices.to.fix = no.of.last.indices.to.fix,
                                                     indices.to.fix = indices.to.fix, lags = lags, extreg = extreg, n.best.extreg	= n.best.extreg,
                                                     use.data.as.ext = use.data.as.ext, lag.externals = lag.externals, consider.as.missing = consider.as.missing,
                                                     whole.period.missing.only = whole.period.missing.only, min.val = min.val, max.val = max.val, digits = digits,
                                                     ...)

  results <- tsrobprep::impute_modelled_data(object = modelled_missings)
}

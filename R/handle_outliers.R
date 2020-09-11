#' Detect and model unreliable outliers of time series data
#'
#' Returns an object of class "tsrobprep" which contains the original data and
#' the modelled values to be imputed. The function handle_outliers detects
#' unreliable outliers given specific threshold and suggests replacement using
#' \link[tsrobprep]{model_missing_data} function. The function analyses both
#' absolute level of the data and the differential series. The function is
#' designed for numerical data only.
#' @param data an input vector, matrix or data frame of dimension nobs x nvars
#' possibly containing unreliable outliers; each column is a variable.
#' @param tau the quantile(s) of the replaced values to be estimated in the
#' quantile regression. Tau accepts all values in (0,1), the default is 0.5.
#' @param S a number of observations per period, e.g. per day.
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
#' @param extreg a vector, matrix or data frame of data containing external
#' regressors.
#' @param use.data.as.ext logical specifying whether to use the remaining
#' variables in the data as external regressors or not.
#' @param min.val a single value or a vector of length nvars providing the
#' minimum possible value of each variable in the data. If a single value, then
#' it applies to all variables. By default set to -Inf.
#' @param max.val a single value or a vector of length nvars providing the
#' maximum possible value of each variable in the data. If a single value, then
#' it applies to all variables. By default set to Inf.
#' @param ... additional arguments for the \link[tsrobprep]{model_missing_data}
#' function.
#' @details The function recognizes two types of outliers: on the absolute level
#' and in the differential series. As outliers are considered all observations
#' that in any of the above mentioned series fall above (below)
#'
#' interquantile_range * margin
#'
#' where margin is given by the user and
#'
#' interquantile_range = upper_quantile - lower_quantile.
#'
#' The upper and lower quantiles are derived from the data based on the
#' outlier.upper.cap and outlier.lower.cap arguments, respectively. In the case
#' of outliers in the differential series, the function uses also smoothing, i.e. if
#' after replacing of the outliers in the differential series, the series still
#' exhibits outliers, the neighbouring observations are also replaced in order
#' to smooth the "jump". The outlier values are replaced using the
#' \link[tsrobprep]{model_missing_data} function.
#'
#' The modelled values can be imputed using \link[tsrobprep]{impute_modelled_data}
#' function.
#' @return An object of class "tsrobprep" which contains the original data, the
#' indices of the data that were modelled, the given quantile values, a list of
#' sparse matrices with the modelled data to be imputed and a list of the numbers
#' of models estimated for every variable.
#' @examples
#' outliers.handled <- handle_outliers(data = GBload[,-1], S = 48, tau = 0.5,
#'                     no.of.last.indices.to.fix = dim(GBload)[1], min.val = 0)
#' outliers.handled$estimated.models
#' outliers.handled$replaced.indices
#' new.GBload <- impute_modelled_data(outliers.handled)
#' @export
#' @seealso \link[tsrobprep]{model_missing_data}, \link[tsrobprep]{impute_modelled_data},
#' \link[tsrobprep]{auto_data_cleaning}
handle_outliers <- function(data, tau = 0.5, S, no.of.last.indices.to.fix = S, indices.to.fix = NULL, outlier.lower.cap = 2*log(dim(as.matrix(data))[1]),
                            outlier.upper.cap = log(dim(as.matrix(data))[1]), margin = outlier.upper.cap, max.periods.to.smooth = 48, extreg = NULL, use.data.as.ext = FALSE,
                            min.val = -Inf, max.val = Inf, ...){
  data.original <- as.matrix(data)
  ###validate the variables' correctness - basic validation
  if(dim(data.original)[1]==0) stop("provided data is of length 0") else nobs <- dim(data.original)[1]; nvars <- dim(data.original)[2]
  if(!is.numeric(tau) | min(tau) <=0 | max(tau) >= 1) stop("provided tau is incorrect")
  if(!(S == floor(S) & S>0)) stop("provided S is incorrect")
  if(no.of.last.indices.to.fix < 1 | no.of.last.indices.to.fix > nobs) stop("provided no.of.last.indices.to.fix is incorrect")
  if(!is.null(indices.to.fix) & !all(is.element(indices.to.fix, 1:nobs))) stop("provided indices exceed the index set")
  if(!is.numeric(outlier.lower.cap)) stop("provided outlier.lower.cap is incorrect")
  if(!is.numeric(outlier.upper.cap)) stop("provided outlier.upper.cap is incorrect")
  if(!is.numeric(margin)) stop("provided margin is incorrect")
  if(!is.numeric(max.periods.to.smooth)) stop("provided max.periods.to.smooth is incorrect")
  if(!is.null(extreg) & (is.numeric(extreg) | is.matrix(extreg) | is.data.frame(extreg))) if(nobs != dim(as.matrix(extreg))[1])
    stop("provided data and external regressors are of different number of observations")
  if(!is.logical(use.data.as.ext)) stop("provided use.data.as.ext is not logical")
  if(length(min.val)!=1 & length(min.val) != nvars) stop("provided min.val is not correct")
  if(length(max.val)!=1 & length(max.val) != nvars) stop("provided max.val is not correct")


  #define the indices to be fixed
  if(is.null(indices.to.fix)) indices.to.fix <- 1:no.of.last.indices.to.fix + nobs - no.of.last.indices.to.fix

  ## define the min/max values
  if(length(min.val)==1) min.vals <- rep(min.val, nvars) else min.vals <- min.val
  if(length(max.val)==1) max.vals <- rep(max.val, nvars) else max.vals <- max.val

  extreg.original <- extreg
  handled.outliers <- list()
  replaced.values <- list()
  models.no <- list()

  for(var.id in 1:nvars){
    handled.outliers[[var.id]] <- list()
    data.to.model <- data.original[,var.id]
    min.val <- min.vals[var.id]
    max.val <- max.vals[var.id]
    if(use.data.as.ext & nvars>1){
      extreg.to.model <- cbind(extreg.original, data.original[,-var.id])
    } else extreg.to.model <- extreg.original

    # Step 1: seek for unreliable outliers
    upper.quants <- stats::quantile(data.to.model, probs = c((nobs - outlier.upper.cap)/nobs, 1 - (nobs - outlier.upper.cap)/nobs), na.rm = TRUE)
    lower.quants <- stats::quantile(data.to.model, probs = c((nobs - outlier.lower.cap)/nobs, 1 - (nobs - outlier.lower.cap)/nobs), na.rm = TRUE)
    interquant.range <- upper.quants - lower.quants
    med <- stats::median(data.to.model, na.rm = TRUE)
    outliers <- integer()
    if(upper.quants[1] + interquant.range[1]*margin == med) warning("the top margin equals the median of the data") else outliers <-
      intersect(which(data.to.model > upper.quants[1] + interquant.range[1]*margin), indices.to.fix)
    if(upper.quants[2] + interquant.range[2]*margin == med) warning("the bottom margin equals the median of the data") else outliers <-
      c(outliers, intersect(which(data.to.model < upper.quants[2] + interquant.range[2]*margin), indices.to.fix))
    handled.outliers[[var.id]][["absolute"]] <- outliers
    data.to.model[outliers] <- NA


    # Step 2: seek for unreliable jumps
    jump.outliers <- numeric()
    k <- 0
    iter <- 0
    while(TRUE){
      jumps <- diff(data.to.model)
      n <- length(jumps)
      upper.quants <- stats::quantile(jumps, probs = c((n - outlier.upper.cap)/n, 1 - (n - outlier.upper.cap)/n), na.rm = T)
      lower.quants <- stats::quantile(jumps, probs = c((n - outlier.lower.cap)/n, 1 - (n - outlier.lower.cap)/n), na.rm = T)
      interquant.range <- upper.quants - lower.quants
      med <- stats::median(jumps, na.rm = T)
      outliers.from.jumps <- numeric()
      if(upper.quants[1] + interquant.range[1]*margin == med) warning("the top margin equals the median of the differences") else
        outliers.from.jumps <- intersect(which(jumps > upper.quants[1] + interquant.range[1]*margin), indices.to.fix)
      if(upper.quants[2] + interquant.range[2]*margin == med) warning("the bottom margin equals the median of the differences") else
        outliers.from.jumps <- c(outliers.from.jumps, intersect(which(jumps < upper.quants[2] + interquant.range[2]*margin), indices.to.fix))
      if(length(outliers.from.jumps)==0) break
      outliers.temp <- unique(c(outliers.from.jumps, outliers.from.jumps+1))
      outliers.test <- outliers.temp %in% jump.outliers
      while(any(outliers.test)){
        new.outliers <- unique(setdiff(as.numeric(sapply(outliers.temp[outliers.test], function(x) x + -2^(iter-1):2^(iter-1))), outliers.temp))
        outliers.temp <- unique(c(outliers.temp, new.outliers))
        outliers.test <- new.outliers %in% jump.outliers
      }
      new.outliers <- setdiff(outliers.temp, jump.outliers)
      jump.outliers[1:length(new.outliers) + k] <- new.outliers
      data.to.model[outliers.temp] <- NA
      new.data <- tsrobprep::model_missing_data(data = data.to.model, tau = 0.5, S = S, use.data.as.ext =FALSE, indices.to.fix = sort(outliers.temp),
                                                extreg = extreg.to.model, min.val = min.val, max.val = max.val, ...)
      data.to.model <- tsrobprep::impute_modelled_data(object = new.data, tau = 0.5)[,1]
      if(iter>log(max.periods.to.smooth, 2)){
        warning(paste0("Variable ", var.id, ": Reached maximum number of days to smooth. Please check if the data is trustworthy."))
        break
        }
      k <- k + length(new.outliers)
      iter <- iter + 1
    }
    handled.outliers[[var.id]][["difference"]] <- sort(jump.outliers)

    ind.to.fix <- sort(c(outliers, jump.outliers))
    data.to.model[ind.to.fix] <- NA
    new.data <- tsrobprep::model_missing_data(data = data.to.model, tau = tau, S = S, use.data.as.ext =FALSE, indices.to.fix = ind.to.fix,
                                              extreg = extreg.to.model, min.val = min.val, max.val = max.val, ...)
    replaced.values[[var.id]] <- new.data$replaced.values[[1]]
    models.no[[var.id]] <- new.data$estimated.models[[1]]
  }
  result <- list(original.data = data.original, replaced.indices =  handled.outliers, tau = tau,
                 replaced.values = replaced.values, estimated.models = models.no)
  class(result) <- "tsrobprep"
  return(result)
}


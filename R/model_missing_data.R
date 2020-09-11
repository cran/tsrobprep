#' Model missing time series data
#'
#' Returns an object of class "tsrobprep" which contains the original data and
#' the modelled missing values to be imputed. The function model_missing_data
#' models missing values in a time series data using the quantile regression
#' implemented in quantreg package. The model uses autoregression on the time
#' series as explanatory variables as well as the provided external variables.
#' The function is designed for numerical data only.
#' @param data an input vector, matrix or data frame of dimension nobs x nvars
#' containing missing values; each column is a variable.
#' @param tau the quantile(s) of the missing values to be estimated in the
#' quantile regression. Tau accepts all values in (0,1), the default is 0.5.
#' @param S a number of observations per period, e.g. per day.
#' @param no.of.last.indices.to.fix a number of observations in the tail of
#' the data to be fixed, by default set to S.
#' @param indices.to.fix indices of the data to be fixed. If NULL, then it is
#' calculated based on the no.of.last.indices.to.fix parameter. Otherwise, the
#' no.of.last.indices.to.fix parameter is ignored.
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
#' @details The function uses quantile regression in order to model missing
#' values and prepare it for imputation. In this purpose the
#' \link[quantreg]{rq.fit.fnb} function from quantreg package is used. The
#' function computes the quantile regression methods utilizing the Frisch-Newton
#' algorithm for user-specified quantile values. The modelled values can be
#' imputed using \link[tsrobprep]{impute_modelled_data} function.
#' @return An object of class "tsrobprep" which contains the original data, the
#' indices of the data that were modelled, the given quantile values, a list of
#' sparse matrices with the modelled data to be imputed and a list of the numbers
#' of models estimated for every variable.
#' @examples
#' model.miss <- model_missing_data(data = GBload[,-1], S = 48, tau = 0.5,
#'               no.of.last.indices.to.fix = dim(GBload)[1], consider.as.missing = 0,
#'               min.val = 0)
#' model.miss$estimated.models
#' model.miss$replaced.indices
#' new.GBload <- impute_modelled_data(model.miss)
#' @export
#' @seealso \link[tsrobprep]{impute_modelled_data}, \link[tsrobprep]{handle_outliers},
#' \link[tsrobprep]{auto_data_cleaning}
model_missing_data <- function(data, tau = 0.5, S, no.of.last.indices.to.fix = S, indices.to.fix = NULL,
                               lags = c(1,2,S,S+1, 7*S, 7*S+1, -1, -2, -S, -S-1, -7*S, -7*S-1),
                               extreg = NULL, n.best.extreg = NULL, use.data.as.ext = FALSE,
                               lag.externals = FALSE, consider.as.missing = NULL,
                               whole.period.missing.only = FALSE, min.val = -Inf, max.val = Inf,
                               digits = 3, ...){

  data.original <- as.matrix(data)
  ###validate the variables' correctness - basic validation
  if(dim(data.original)[1]==0) stop("provided data is of length 0") else nobs <- dim(data.original)[1]; nvars <- dim(data.original)[2]
  if(!is.numeric(tau) | min(tau) <=0 | max(tau) >= 1) stop("provided tau is incorrect")
  if(!(S == floor(S) & S>0)) stop("provided S is incorrect")
  if(no.of.last.indices.to.fix < 1 | no.of.last.indices.to.fix > nobs) stop("provided no.of.last.indices.to.fix is incorrect")
  if(!is.null(indices.to.fix) & !all(is.element(indices.to.fix, 1:nobs))) stop("provided indices exceed the index set")
  if(!is.numeric(lags)) stop("provided lags is not numeric")
  if(!is.null(extreg) & (is.numeric(extreg) | is.matrix(extreg) | is.data.frame(extreg))) if(nobs != dim(as.matrix(extreg))[1])
    stop("provided data and external regressors are of different number of observations")
  if(!is.null(n.best.extreg)) if(!(n.best.extreg == floor(n.best.extreg) & n.best.extreg>0)) stop("provided n.best.extreg is incorrect")
  if(!is.logical(use.data.as.ext)) stop("provided use.data.as.ext is not logical")
  if(!is.logical(lag.externals)) stop("provided lag.externals is not logical")
  if(!is.logical(whole.period.missing.only)) stop("provided whole.period.missing.only is not logical")
  if(length(min.val)!=1 & length(min.val) != nvars) stop("provided min.val is not correct")
  if(length(max.val)!=1 & length(max.val) != nvars) stop("provided max.val is not correct")

  #define the indices to be fixed
  if(is.null(indices.to.fix)) indices.to.fix <- 1:no.of.last.indices.to.fix + nobs - no.of.last.indices.to.fix

  ## define the min/max values
  if(length(min.val)==1) min.vals <- rep(min.val, nvars) else min.vals <- min.val
  if(length(max.val)==1) max.vals <- rep(max.val, nvars) else max.vals <- max.val


  extreg.original <- extreg
  considered.missing <- list()
  replaced.values <- list()
  models.no <- list()
  for(var.id in 1:nvars){
    considered.missing[[var.id]] <- list()
    data.to.model <- data.original[,var.id]
    min.val <- min.vals[var.id]
    max.val <- max.vals[var.id]
    #join the extreg and remaining variables of data
    if(use.data.as.ext & nvars>1){
      extreg.to.model <- cbind(extreg.original, data.original[,-var.id])
    } else extreg.to.model <- extreg.original

    #remove values considered as missings from externals
    if(!is.null(extreg.to.model)){
      extreg.to.model <- as.matrix(extreg.to.model)
      if(!is.null(consider.as.missing) & !whole.period.missing.only){
        for(i in 1:dim(extreg.to.model)[2]){
          ZEROs <- rep(FALSE, nobs)
          ZEROs[indices.to.fix] <- extreg.to.model[indices.to.fix,i] %in% consider.as.missing
          ZEROs[is.na(ZEROs)] <- FALSE
          extreg.to.model[ZEROs,i] <- NA
        }
      }
      if(!is.null(consider.as.missing) & whole.period.missing.only){
        for(i in 1:dim(extreg.to.model)[2]){
          day.by.day.matrix <- matrix(extreg.to.model[,i] %in% consider.as.missing, nrow = S, byrow = F)
          whole.period.missings <- which(apply(day.by.day.matrix, 2, all))
          whole.period.indices <- which(col(day.by.day.matrix) %in% whole.period.missings)
          ZEROs <- rep(FALSE, nobs)
          ZEROs[intersect(whole.period.indices, indices.to.fix)] <- TRUE
          extreg.to.model[ZEROs,i] <- NA
        }
      }
    }




    #if threshold specified, calculate the correlation level between the data and external regressors and use only the ones that exceed the threshold
    if(!is.null(n.best.extreg) & !is.null(extreg.to.model)){
      cor.ext <- order(abs(as.numeric(stats::cor(data.to.model, extreg.to.model, use = "pairwise.complete.obs"))), decreasing = T)
      extreg.to.model <- as.matrix(extreg.to.model)[,utils::head(cor.ext, n.best.extreg)]
    }

    #function for lagging the data
    get.lagged<- function(lag, Z){
      if(lag>=0) c( rep(NA, lag),  Z[(1+lag):(length(Z)) -lag ]  ) else c(Z[-(1:(-lag))], rep(NA, -lag) )
    }

    lags <- unique(lags)

    ## lag externals
    if(lag.externals){
      new.extreg <- extreg.to.model
      for(ext.col in 1:dim(extreg.to.model)[2]){
        new.extreg <- cbind(new.extreg, sapply(lags, get.lagged, Z = extreg.to.model[,ext.col]))
      }
      extreg.to.model <- new.extreg
    }


    ### create a vector to iterate over
    #indices to be fixed cause NA
    NAs <- rep(FALSE, nobs)
    NAs[indices.to.fix] <- is.na(data.to.model[indices.to.fix])
    considered.missing[[var.id]][["NA"]] <- which(NAs)
    #indices to be fixed cause given values
    if(!is.null(consider.as.missing) & !whole.period.missing.only){
      for(i in consider.as.missing){
        ZEROs <- rep(FALSE, nobs)
        ZEROs[indices.to.fix] <- data.to.model[indices.to.fix] %in% i
        ZEROs[is.na(ZEROs)] <- FALSE
        considered.missing[[var.id]][[as.character(i)]] <- which(ZEROs)
        data.to.model[considered.missing[[var.id]][[as.character(i)]]] <- NA
        NAs <- NAs | ZEROs
      }
    }
    if(!is.null(consider.as.missing) & whole.period.missing.only){
      for(i in consider.as.missing){
        day.by.day.matrix <- matrix(data.to.model %in% i, nrow = S, byrow = F)
        whole.period.missings <- which(apply(day.by.day.matrix, 2, all))
        whole.period.indices <- which(col(day.by.day.matrix) %in% whole.period.missings)
        ZEROs <- rep(FALSE, nobs)
        ZEROs[intersect(whole.period.indices, indices.to.fix)] <- TRUE
        considered.missing[[var.id]][[as.character(i)]] <- which(ZEROs)
        data.to.model[considered.missing[[var.id]][[as.character(i)]]] <- NA
        NAs <- NAs | ZEROs
      }
    }


    #lag the data
    X <- sapply(lags, get.lagged, Z = data.to.model)
    #bind the lagged and external regressors
    if(!is.null(extreg.to.model)) X <- cbind(X, extreg.to.model)
    df.orig <- cbind(y = data.to.model, X, intercept = 1)
    tau.len <- length(tau)
    df.new <- Matrix::Matrix(0, nrow = nobs, ncol = tau.len)
    df.learn <- df.orig[which( !is.na(df.orig[,1])) ,]
    bad.regressors <- apply(is.na(df.learn), 2, mean) < 0.5
    df.learn <- df.learn[, bad.regressors]
    df.orig <- df.orig[,bad.regressors]
    df.learn <- df.learn[which( !apply(is.na(df.learn),1, any) ),]
    # p is the number of lags
    p <- length(lags)


    #if no NAs, no procedure
    if(any(NAs)){
      iter <- which(NAs)

      # iterating over the missing points and interpolating them

      for(tau.no in 1:tau.len){
        models <- list()
        cols.to.remove <- list()
        for(i.miss in 1:length(iter)){
          missing.index <- iter[i.miss]
          #cat(i.miss/length(iter), "\n")
          test.data <- df.orig[missing.index,]
          test <- !is.na(test.data)
          test.str <- paste(as.numeric(test), collapse="")
          if(is.null(models[[test.str]])){
            models[[test.str]] <- try({quantreg::rq.fit.fnb(x = as.matrix(df.learn[,test]), y = df.learn[,1], tau = tau[tau.no], ...)$coefficients},
                                      silent = TRUE)
            while(class(models[[test.str]]) == "try-error"){
              error_type <- attr(models[[test.str]],"condition")
              if(error_type$message!= "Error info =  1 in stepy: singular design") stop(error_type)
              warning("singular design in regression matrix")
              if(length(cols.to.remove[[test.str]])==0){
                cols.to.remove[[test.str]] <- c(cols.to.remove[[test.str]],
                                                which.max(rowSums((abs(svd(df.learn[,test])$v) < .Machine$double.eps) + 0)))
              } else cols.to.remove[[test.str]] <- c(cols.to.remove[[test.str]],
                                                     which.max(rowSums((abs(svd(df.learn[,test][, - cumsum(cols.to.remove[[test.str]])])$v) < .Machine$double.eps) + 0)))

              models[[test.str]] <- try({quantreg::rq.fit.fnb(x = as.matrix(df.learn[,test][, -cumsum(cols.to.remove[[test.str]])]),
                                                                   y = df.learn[,1], tau = tau[tau.no], ...)$coefficients}, silent = TRUE)
            }
            #models[[test.str]] <- quantreg::rq.fit.fnb(x = as.matrix(df.learn[,test]), y = df.learn[,1], tau = tau[tau.no], ...)$coefficients
            #models[[test.str]] <- quantreg::rq.fit.fnb(x = as.matrix(df.learn[,test]), y = df.learn[,1], tau = tau[tau.no])$coefficients
          }
          #if(is.null(models[[test.str]])) models[[test.str]] <- quantreg::rq.fit.fnb(x = as.matrix(df.orig[,test]), y = df.learn[,1], tau = tau[tau.no])$coefficients
          if(length(cols.to.remove[[test.str]])==0){
            val <- round(pmin(pmax(as.numeric(test.data[test] %*% models[[test.str]]),min.val),max.val), digits)
          } else
            val <- round(pmin(pmax(as.numeric(test.data[test][- cumsum(cols.to.remove[[test.str]])] %*% models[[test.str]]),min.val),max.val), digits)


          df.new[missing.index,tau.no] <- val

        }
      }

      models.no[[var.id]] <- length(models)
      replaced.values[[var.id]] <- df.new
    }
  }
  result <- list(original.data = data.original, replaced.indices =  considered.missing, tau = tau,
                 replaced.values = replaced.values, estimated.models = models.no)
  class(result) <- "tsrobprep"
  return(result)
}
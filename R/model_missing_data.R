#' Model missing time series data
#'
#' Returns an object of class "tsrobprep" which contains the original data and
#' the modelled missing values to be imputed. The function model_missing_data
#' models missing values in a time series data using a robust time series
#' decomposition with the weighted lasso or the quantile regression. The model
#' uses autoregression on the time series as explanatory variables as well as
#' the provided external variables. The function is designed for numerical data
#' only.
#' @param data an input vector, matrix or data frame of dimension nobs x nvars
#' containing missing values; each column is a variable.
#' @param S a number or vector describing the seasonalities (S_1, ..., S_K) in
#' the data, e.g. c(24, 168) if the data consists of 24 observations per day
#' and there is a weekly seasonality in the data.
#' @param tau the quantile(s) of the missing values to be estimated in the
#' quantile regression. Tau accepts all values in (0,1). If NULL, then the
#' weighted lasso regression is performed.
#' @param no.of.last.indices.to.fix a number of observations in the tail of
#' the data to be fixed, by default set to first element of S.
#' @param indices.to.fix indices of the data to be fixed. If NULL, then it is
#' calculated based on the no.of.last.indices.to.fix parameter. Otherwise, the
#' no.of.last.indices.to.fix parameter is ignored.
#' @param replace.recursively if TRUE then the algorithm uses replaced values
#' to model the remaining missings.
#' @param p a number or vector of length(S) = K indicating the order of a
#' K-seasonal autoregressive process to be estimated. If NULL, chosen
#' data-based.
#' @param mirror if TRUE then autoregressive lags up to order p are not only
#' added to the seasonalities but also subtracted.
#' @param lags a numeric vector with the lags to use in the autoregression.
#' Negative values are accepted and then also the "future" observations are
#' used for modelling. If not NULL, p and mirror are ignored.
#' @param extreg a vector, matrix or data frame of data containing external
#' regressors; each column is a variable.
#' @param n.best.extreg a numeric value specifying the maximal number of
#' considered best correlated external regressors (selected in decreasing
#' order). If NULL, then all variables in extreg are used for modelling.
#' @param use.data.as.ext logical specifying whether to use the remaining
#' variables in the data as external regressors or not.
#' @param lag.externals logical specifying whether to lag the external
#' regressors or not. If TRUE, then the algorithm uses the lags specified in
#' parameter lags.
#' @param consider.as.missing a vector of numerical values which are considered
#' as missing in the data.
#' @param whole.period.missing.only if FALSE, then all observations which
#' correspond to the values of consider.as.missing are treated as missings. If
#' TRUE, then only consecutive observations of specified length are considered
#' (length is defined by first element of S).
#' @param debias if TRUE, the recursive replacement is additionally debiased.
#' @param min.val a single value or a vector of length nvars providing the
#' minimum possible value of each variable in the data. If a single value, then
#' it applies to all variables. By default set to -Inf.
#' @param max.val a single value or a vector of length nvars providing the
#' maximum possible value of each variable in the data. If a single value, then
#' it applies to all variables. By default set to Inf.
#' @param Cor_thres a single value providing the correlation threshold from
#' which external regressors are considered in the quantile regression.
#' @param digits integer indicating the number of decimal places allowed
#' in the data, by default set to 3.
#' @param ICpen is the information criterion penalty for lambda choice in the
#' \link[glmnet]{glmnet} algorithm. It can be a string: "BIC", "HQC" or "AIC",
#' or a fixed number.
#' @param decompose.pars named list containing additional arguments for the
#' \link[tsrobprep]{robust_decompose} function.
#' @param ... additional arguments for the \link[glmnet]{glmnet} or
#' \link[quantreg]{rq.fit.fnb} algorithms.
#' @details The function uses robust time series decomposition with weighted
#' lasso or quantile regression in order to model missing values and prepare it
#' for imputation. In this purpose the \link[tsrobprep]{robust_decompose}
#' function together with the \link[glmnet]{glmnet} are used in case of mean
#' regression, i.e. tau = NULL. In case of quantile regression, i.e.
#' tau != NULL the \link[tsrobprep]{robust_decompose} function is used together
#' with the \link[quantreg]{rq.fit.fnb} function. The modelled values can be
#' imputed using \link[tsrobprep]{impute_modelled_data} function.
#' @return An object of class "tsrobprep" which contains the original data, the
#' indices of the data that were modelled, the given quantile values, a list of
#' sparse matrices with the modelled data to be imputed and a list of the
#' numbers of models estimated for every variable.
#' @examples
#' \dontrun{
#' model.miss <- model_missing_data(
#'     data = GBload[,-1], S = c(48,7*48),
#'     no.of.last.indices.to.fix = dim(GBload)[1], consider.as.missing = 0,
#'     min.val = 0
#' )
#' model.miss$estimated.models
#' model.miss$replaced.indices
#' new.GBload <- impute_modelled_data(model.miss)
#' }
#' @export
#' @seealso \link[tsrobprep]{robust_decompose},
#' \link[tsrobprep]{impute_modelled_data}, \link[tsrobprep]{detect_outliers},
#' \link[tsrobprep]{auto_data_cleaning}
model_missing_data <- function(
  data, S, tau = NULL, no.of.last.indices.to.fix = S[1], indices.to.fix = NULL,
  replace.recursively = TRUE, p = NULL, mirror = FALSE, lags = NULL,
  extreg = NULL, n.best.extreg = NULL, use.data.as.ext = FALSE,
  lag.externals = FALSE, consider.as.missing = NULL,
  whole.period.missing.only = FALSE, debias = FALSE, min.val = -Inf,
  max.val = Inf, Cor_thres = 0.5, digits = 3, ICpen = "BIC",
  decompose.pars = list(), ...){

  # save the data as matrix
  data.original <- as.matrix(data)
  # validate the variables' correctness - basic validation
  if(dim(data.original)[1]==0) stop("provided data is of length 0") else
    nobs <- dim(data.original)[1]; nvars <- dim(data.original)[2]
  if(!is.null(tau)){
    if(min(tau) <=0 | max(tau) >= 1) stop("provided tau is incorrect")
  } else tau <- "mean"
  if(no.of.last.indices.to.fix < 1 | no.of.last.indices.to.fix > nobs)
    stop("provided no.of.last.indices.to.fix is incorrect")
  if(!is.null(indices.to.fix) & !all(is.element(indices.to.fix, 1:nobs)))
    stop("provided indices are not part of the index set")
  if(!is.null(extreg) & (is.numeric(extreg) | is.matrix(extreg) |
                         is.data.frame(extreg)))
    if(nobs != dim(as.matrix(extreg))[1])
      stop(paste0("provided data and external regressors are of different",
                  " number of observations"))
  if(!is.null(n.best.extreg))
    if(!(n.best.extreg == floor(n.best.extreg) & n.best.extreg>0))
      stop("provided n.best.extreg is incorrect")
  if(!is.logical(use.data.as.ext))
    stop("provided use.data.as.ext is not logical")
  if(!is.logical(lag.externals))
    stop("provided lag.externals is not logical")
  if(!is.logical(whole.period.missing.only))
    stop("provided whole.period.missing.only is not logical")
  if(length(min.val)!=1 & length(min.val) != nvars)
    stop("provided min.val is not correct")
  if(length(max.val)!=1 & length(max.val) != nvars)
    stop("provided max.val is not correct")

  #define the indices to be fixed
  if(is.null(indices.to.fix)) indices.to.fix <-
      seq_len(no.of.last.indices.to.fix) + nobs - no.of.last.indices.to.fix

  ## define the min/max values
  if(length(min.val)==1) min.vals <- rep(min.val, nvars) else
    min.vals <- min.val
  if(length(max.val)==1) max.vals <- rep(max.val, nvars) else
    max.vals <- max.val

  # save orignal extreg
  extreg.original <- extreg

  # define output lists
  considered.missing <- list()
  replaced.values <- list()
  models.no <- list()

  # we work in a loop over variables (nvars)
  for(var.id in seq_len(nvars)){
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

      # the case of whole.period.missing.only = FALSE
      if(!is.null(consider.as.missing) & !whole.period.missing.only){
        for(i in 1:dim(extreg.to.model)[2]){
          ZEROs <- rep(FALSE, nobs)
          ZEROs[indices.to.fix] <-
            extreg.to.model[indices.to.fix,i] %in% consider.as.missing
          ZEROs[is.na(ZEROs)] <- FALSE
          extreg.to.model[ZEROs,i] <- NA
        }
      }
      # the case of whole.period.missing.only = TRUE
      if(!is.null(consider.as.missing) & whole.period.missing.only){
        for(i in 1:dim(extreg.to.model)[2]){
          day.by.day.matrix <- matrix(extreg.to.model[,i] %in%
                                        consider.as.missing, nrow = S[1],
                                      byrow = F)
          whole.period.missings <- which(apply(day.by.day.matrix, 2, all))
          whole.period.indices <-
            which(col(day.by.day.matrix) %in% whole.period.missings)
          ZEROs <- rep(FALSE, nobs)
          ZEROs[intersect(whole.period.indices, indices.to.fix)] <- TRUE
          extreg.to.model[ZEROs,i] <- NA
        }
      }
    }

    # function preparing the lags
    get.lags <- function(S=24, p=1, mirror=FALSE, lags=NULL){
      if(is.null(lags)){
        SS <- unique(c(1,S))
        K <- length(SS)

        p <- rep(p, length.out=K) ## recycle p
        Slags <- list()
        for(i.k in 1:K){
          Slags[[i.k]] <- list()
          for(i.p in 1:(p[i.k]+1)) {
            tmp <- i.p-1
            if(mirror) tmp <- c(-tmp,tmp)
            Slags[[i.k]][[i.p]] <- SS[i.k] +tmp
          }
        }
        # exp_lags <- c()
        # for(i.k in 1:(K-1)){
        #   tmp <- 2^(1:floor(log(SS[i.k+1]/SS[i.k],base=2)))*SS[i.k]
        #   exp_lags <- c(exp_lags,tmp)
        # }
        #
        # Srec <- c(unlist(Slags), exp_lags)
        Srec <- unlist(Slags)
        Srec <- sort(unique(Srec[Srec>0]))
        lags <- c(-rev(Srec),Srec)
      }
      lags
    }

    # if p not specified, calculate it based on the available data length
    p_new <- ifelse(is.null(p), pmax(
      floor(log(sum(!is.na(data.to.model[indices.to.fix])),10)-1),1), p)

    # derive the lags to be used in the regression
    lags_new <- get.lags(S=S,p=p_new, mirror=mirror, lags = lags)

    #function for lagging the data
    get.lagged <- function(lag, Z){
      if(lag>=0) c(rep(NA, lag),  Z[(1+lag):(length(Z)) -lag ]) else
        c(Z[-(1:(-lag))], rep(NA, -lag) )
    }

    # no duplicates in lags
    lags_new <- unique(lags_new)

    # if max lags exceeding nobs, throw a warning
    if(max(abs(lags_new))>=nobs){
      warning(
        paste0("Some of the given lags were exceeding the number of",
               " observations. The set of lags was accordingly truncated.")
      )
      lags_new <- lags_new[abs(lags_new) < nobs]
    }
    lag_len <- length(lags_new)

    # if lag.externals = TRUE apply the below
    if(lag.externals){
      new.extreg <- extreg.to.model
      for(ext.col in 1:dim(extreg.to.model)[2]){
        new.extreg <- cbind(new.extreg, sapply(lags_new, get.lagged,
                                               Z = extreg.to.model[,ext.col]))
      }
      extreg.to.model <- new.extreg
    }

    # if threshold specified, calculate the correlation level between the data
    # and external regressors and use only the ones that exceed the threshold
    if(!is.null(n.best.extreg) & !is.null(extreg.to.model)){
      cor.ext <- abs(as.numeric(stats::cor(
        data.to.model, extreg.to.model, use = "pairwise.complete.obs")))
      cor.ext_thres <- which(cor.ext >= Cor_thres)
      cor.ext_order <- order(cor.ext, decreasing = T)
      extreg.to.model <- as.matrix(extreg.to.model)[,intersect(
        utils::head(cor.ext_order, n.best.extreg),cor.ext_thres)]
    }




    # create a vector to iterate over
    # indices to be fixed cause NA
    NAs <- rep(FALSE, nobs)
    NAs[indices.to.fix] <- is.na(data.to.model[indices.to.fix])
    considered.missing[[var.id]][["NA"]] <- which(NAs)

    # indices to be fixed cause given consider as missing values
    # the case of whole.period.missing.only = FALSE
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
    # the case of whole.period.missing.only = TRUE
    if(!is.null(consider.as.missing) & whole.period.missing.only){
      for(i in consider.as.missing){
        day.by.day.matrix <- matrix(data.to.model %in% i, nrow = S[1],
                                    byrow = F)
        whole.period.missings <- which(apply(day.by.day.matrix, 2, all))
        whole.period.indices <- which(col(day.by.day.matrix) %in%
                                        whole.period.missings)
        ZEROs <- rep(FALSE, nobs)
        ZEROs[intersect(whole.period.indices, indices.to.fix)] <- TRUE
        considered.missing[[var.id]][[as.character(i)]] <- which(ZEROs)
        data.to.model[considered.missing[[var.id]][[as.character(i)]]] <- NA
        NAs <- NAs | ZEROs
      }
    }

    # run the robust decomposition
    data_decomposed <- do.call(tsrobprep::robust_decompose, c(
      list(x = data.to.model, S = S, extreg = extreg.to.model),
      decompose.pars))

    # remove constant components
    data_decomposed$components <-
      data_decomposed$components[,apply(data_decomposed$components, 2,
                                        stats::var, na.rm=TRUE) != 0, drop = F]

    # subtract the fit of the decomposition
    data.to.model <- data.to.model - data_decomposed$fit

    #lag the data
    X <- sapply(lags_new, get.lagged, Z = data.to.model)
    colnames(X) <- paste0("lag_", lags_new)

    # #bind the lagged and external regressors
    # if(!is.null(extreg.to.model)) X <- cbind(extreg.to.model, X)

    # cbind all variables
    df.orig <- cbind(y = data.to.model, intercept = 1)
    if(dim(data_decomposed$components)[2]>0){
      df.orig <- cbind(df.orig, data_decomposed$components)
    }
    df.orig <- cbind(df.orig, X)

    # get the number of taus
    tau.len <- length(tau)

    # define the output for replaced values
    df.new <- Matrix::Matrix(0, nrow = nobs, ncol = tau.len)

    # get the non-na observations
    act_ind <- !is.na(df.orig[,1])

    # define the learning data frame
    df.learn <- df.orig[act_ind, ]

    # define weights
    df_weights_orig <- data_decomposed$weights[act_ind]

    ### TODO Change static threshold to dynamic threshold
    # if any of the regressors is NA more often than in thr*100% cases then we
    # exclude it from the df.learn
    thr <- 0.5
    bad.regressors <- apply(is.na(df.learn), 2, mean) <= thr
    df.learn <- df.learn[, bad.regressors]
    act_ind <- !apply(is.na(df.learn),1, any)
    df.learn <- df.learn[act_ind ,, drop = FALSE]
    df_weights <- df_weights_orig[act_ind]

    # if too little observations, lower the threshold
    while(nrow(df.learn)<ncol(df.learn)){
      df.learn <- df.orig[which( !is.na(df.orig[,1])) ,]
      thr <- pmax(thr - 0.05,0)
      bad.regressors <- apply(is.na(df.learn), 2, mean) <= thr
      df.learn <- df.learn[, bad.regressors]
      act_ind <- !apply(is.na(df.learn),1, any)
      df.learn <- df.learn[act_ind ,, drop = FALSE]
      df_weights <- df_weights_orig[act_ind]
    }
    df.orig <- df.orig[,bad.regressors]
    reg_num <- dim(df.orig)[2]
    # TODO better
    if(thr < 0.5) warning(paste0("The regression matrix consists of many NAs.",
                                 " The replacement may be inaccurate."))

    #if no NAs, no procedure
    if(any(NAs)){
      if(replace.recursively){
        #create vector to iterate over in such a manner that in bigger gaps
        # the algorithm goes alternately "from left" and "from right"
        iter <- numeric(0)
        gaps <- numeric(0)
        max.counter <- sum(NAs)
        counter <- which.min(NAs) - 1
        na.in.tail <- which.min(NAs[nobs:1]) - 2
        gap_number <- 1
        # if NA in tail, go simply from left to right
        if(na.in.tail>0){
          iter[max.counter + -na.in.tail:0] <- nobs + -na.in.tail:0
          gaps[max.counter + -na.in.tail:0] <-
            gap_number + seq_along(nobs + -na.in.tail:0) - 1
          gap_number <- gap_number + length(nobs + -na.in.tail:0)
          max.counter <- max.counter - na.in.tail - 1
        }
        # if NA in front, go from right to the left
        if(counter>0){
          iter[0:counter] <- counter:1
          new.NAs <- NAs[-(1:counter)]
          gaps[0:counter] <- gap_number + seq_along(counter:1) - 1
          gap_number <- gap_number + length(counter:1)
        } else new.NAs <- NAs
        iter.len <- counter
        # now if any longer gaps appear in the middle of the data set, model
        # alternately "from left" and "from right"
        while(iter.len < max.counter){
          new.start <- which.max(new.NAs)-1
          counter <- counter + new.start
          new.NAs <- new.NAs[-(1:new.start)]
          items.num <- which.min(new.NAs)-1
          # if less than 3 elements in the gap, no point of doing it
          if(items.num < 3){
            iter[(iter.len+1):(iter.len+items.num)] <-
              (counter+1):(counter+items.num)
            gaps[(iter.len+1):(iter.len+items.num)] <-
              gap_number + seq_along((counter+1):(counter+items.num)) - 1
            gap_number <- gap_number + length((counter+1):(counter+items.num))
          } else{
            new.it <- numeric(0)
            forw <- (counter+1):(counter+items.num)
            backw <- (counter+items.num):(counter+1)
            new.it[seq(1, 2*items.num, by = 2)] <- forw
            new.it[seq(2, 2*items.num, by = 2)] <- backw
            iter[(iter.len+1):(iter.len+items.num)] <- new.it[1:items.num]
            gaps[(iter.len+1):(iter.len+items.num)] <- gap_number
            gap_number <- gap_number+1
          }
          new.NAs <- new.NAs[-(1:items.num)]
          counter <- counter + items.num
          iter.len <- iter.len + items.num
        }
        lags.to.replace <- c(0, lags_new)
      } else{
        iter <- which(NAs)
        gaps <- iter
      }

      # derive the IC penalty factor
      if(ICpen=="HQC"){
        ICpen <- 2*log(log(sum(!NAs)))
      } else if (ICpen=="AIC") {
        ICpen <- 2
      } else if (ICpen=="BIC"){
        ICpen <- log(sum(!NAs))
      }

      # iterating over the missing points and interpolating them
      for(tau.no in seq_len(tau.len)){
        # copy of df.orig
        df.orig.copy <- df.orig

        # defining the list of models
        models <- list()
        # and columns to remove in case of singularities
        cols.to.remove <- list()

        for(i.gap in unique(gaps)){
          missing.indices <- sort(iter[gaps == i.gap])

          runs <- ifelse(length(missing.indices)>1, 2, 1)

          val <- matrix(nrow = length(missing.indices), ncol = runs)

          for(i.run in seq_len(runs)){
            indices.temp <- sort(
              missing.indices, decreasing = as.logical((i.run-1)%%2))
            if(runs>1 & debias==TRUE) indices.temp[length(indices.temp)+1] <-
                indices.temp[length(indices.temp)] + (-1)^(i.run-1)
            df.temp <- df.orig.copy

            for(i.miss in seq_along(indices.temp)){
              missing.index <- indices.temp[i.miss]
              # cat(i.miss/length(iter), "\r")

              # choose the data to be used for modelling
              test.data <- df.temp[missing.index,-1]
              if(runs>1 & debias==TRUE){
                lag.ind <- (reg_num + 1 - lag_len):reg_num -1
                test.data[lag.ind[(
                  (i.run-1)*length(lag.ind)/2+1):(i.run*length(lag.ind)/2)]] <-
                  NA
              }
              # check which regressors are available
              test <- !is.na(test.data)
              # convert the test to string
              test.str <- paste(as.numeric(test), collapse="")

              # if no such model was estimated so far, estimate a new one
              # else - use already estimated coefficients
              if(is.null(models[[test.str]])){
                # in the case of tau = NULL, use weighted lasso
                if(tau[1] == "mean"){
                  # if constant, just skip the modelling part
                  if(all(df.learn[,1] == df.learn[1,1])){
                    models[[test.str]] <- "const"
                  } else {
                    # fit the lasso model with n = 100 lambdas
                    model <- glmnet::glmnet(
                      x = as.matrix(df.learn[,-1][,test]), y = df.learn[,1],
                      weights = df_weights, ...)

                    # choose the best model IC-based
                    RSS <- (1 - model$dev.ratio) * model$nulldev
                    IC <- log(RSS) + model$df * ICpen / model$nobs
                    idsel <- which.min(IC)

                    # save the coefficients
                    models[[test.str]] <- model$beta[, idsel]
                    models[[test.str]][1] <- model$a0[idsel]
                  }


                } else { # otherwise, use quantile regression
                  # make warnings count as errors, for a while
                  op <- options(warn=2)
                  # try to fit the model and in case of singular design warning
                  # the model will be estimated using less number of regressors
                  models[[test.str]] <- try({quantreg::rq.fit.fnb(
                    x = as.matrix(df.learn[,-1][,test]), y = df.learn[,1],
                    tau = tau[tau.no])$coefficients}, silent = TRUE)
                  # if error and singular design, delete regressors until no
                  # singular design
                  while(class(models[[test.str]]) == "try-error"){
                    error_type <- attr(models[[test.str]],"condition")
                    if(! grepl(pattern = "singular design",
                               x = error_type$message)) stop(error_type)
                    op <- options(warn=1)
                    # save the columns to be removed
                    if(length(cols.to.remove[[test.str]])==0){
                      cols.to.remove[[test.str]] <- which.max(rowSums((abs(
                        svd(df.learn[,test])$v) < .Machine$double.eps) + 0))
                    } else cols.to.remove[[test.str]] <-
                      c(cols.to.remove[[test.str]],
                        (1:ncol(df.learn[,test]))[ -cols.to.remove[[test.str]]][
                          which.max(rowSums((abs(svd(df.learn[,test][
                          , - cols.to.remove[[test.str]]])$v) <
                            .Machine$double.eps) + 0))])
                    # try again to fit the model
                    op <- options(warn=2)
                    models[[test.str]] <- try({quantreg::rq.fit.fnb(
                      x = as.matrix(df.learn[,-1][,test][
                        , - cols.to.remove[[test.str]]]), y = df.learn[,1],
                      tau = tau[tau.no])$coefficients}, silent = TRUE)
                  }
                  # make warnings count as warnings again
                  op <- options(warn=1)
                }
              }
              # if constant and lasso, use simply the constant val
              if(models[[test.str]][1] == 'const' & tau[1] == 'mean'){
                val[match(missing.index, missing.indices), i.run] <-
                  df.learn[1,1]
              } else{
                # if no singularities, just calculate the replacement value
                if(length(cols.to.remove[[test.str]])==0){
                  val[match(missing.index, missing.indices), i.run] <-
                    as.numeric(test.data[test] %*% models[[test.str]])
                } else{ # otherwise, raise a warning additionally
                  warning(paste0(
                    "To avoid singular matrix design in regression ",
                    "matrix relevant regressors were removed"))
                  val[match(missing.index, missing.indices), i.run] <-
                    as.numeric(test.data[test][
                      - cols.to.remove[[test.str]]] %*%
                        models[[test.str]])
                }
              }

              if(replace.recursively){
                if.in.df <- missing.index+lags.to.replace < dim(df.temp)[1] &
                  missing.index+lags.to.replace > 0
                df.temp[cbind((missing.index+lags.to.replace)[if.in.df],
                              c(1, (reg_num + 1 - lag_len):reg_num)[if.in.df])
                        ] <- val[match(missing.index, missing.indices), i.run]
              }

              if(runs>1 & i.miss == length(indices.temp) & debias==TRUE){
                val.temp <- as.numeric(test.data[test] %*% models[[test.str]])
                true.temp <- df.orig.copy[missing.index,1]
                val[sort(seq_len(length(missing.indices)),
                         decreasing = as.logical((i.run-1)%%2)),i.run] <-
                  val[sort(seq_len(length(missing.indices)),
                           decreasing = as.logical((i.run-1)%%2)),i.run] -
                  seq.int(1, dim(val)[1])/dim(val)[1]*(val.temp-true.temp)
              }
            }
          }

          if(runs == 2){
            val_weight <- seq(1/dim(val)[1], 1 - 1/dim(val)[1],
                              length.out = dim(val)[1])
            val <- val[,1]*(1-val_weight)+val[,2]*val_weight
          }


          df.new[missing.indices, tau.no] <- round(pmin(pmax(val +
              data_decomposed$fit[missing.indices],min.val),max.val), digits)
          # if replacing recursively, one must also replace the values
          # in the regressors
          for(i.miss in seq_along(missing.indices)){
            missing.index <- missing.indices[i.miss]
            if(replace.recursively){
              if.in.df <- missing.index+lags.to.replace < dim(df.new)[1] &
                missing.index+lags.to.replace > 0
              df.orig.copy[cbind(
                (missing.index+lags.to.replace)[if.in.df],
                c(1, (reg_num + 1 - lag_len):reg_num)[if.in.df])
                ] <- as.numeric(val)[i.miss]
            }
          }


        }
      }

      models.no[[var.id]] <- length(models)
      # sort the quantiles in case the values are not increasing
      if(tau.len > 1){
        df.new <- Matrix::Matrix(t(apply(df.new, 1, sort, decreasing = FALSE)))
      }
      replaced.values[[var.id]] <- df.new
    }
  }

  result <- list(
    original.data = data.original, replaced.indices =  considered.missing,
    tau = tau, replaced.values = replaced.values, estimated.models = models.no
  )
  class(result) <- "tsrobprep"
  return(result)
}

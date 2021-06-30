#' Robust time series seasonal decomposition
#'
#' Decompose a time series into trend, level and potentially multiple seasonal
#' components including all interactions. The function allows for missings.
#' @param x a time series.
#' @param S a number or vector describing the seasonalities (S_1, ..., S_K) in
#' the data, e.g. c(24, 168) if the data consists of 24 observations per day
#' and there is a weekly seasonality in the data.
#' @param wsize is filter/rolling med size
#' @param use.trend if TRUE, uses standard decomposition. If FALSE, uses no
#' trend component.
#' @param K a sigma (standard deviation) bound. The observations that exceed
#' sigma*K become reduced weight in the regression.
#' @param ICpen is the information criterion penalty, e.g. string "BIC", "HQC"
#' or "AIC", or a fixed number.
#' @param extreg a vector, matrix or data frame of data containing external
#' regressors; each column is a variable.
#' @param use.autoregressive if TRUE, removes the autoregression from the
#' series. If NULL, it is derived data based.
#' @details
#' @return A list which contains a vector of fitted values, a vector of weights
#' given to the original time series, and a matrix of components of the
#' decomposition.
#' @examples
#' \dontrun{
#' GBload.decomposed <- robust_decompose(GBload[,-1], S = c(48,7*48))
#' head(GBload.decomposed$components)
#' }
#' @export

robust_decompose <- function(
  x, S, wsize= max(2*max(S),25), use.trend=TRUE, K=4, ICpen="BIC", extreg=NULL,
  use.autoregressive = NULL){
# 13 as by forecast::stlf
  # TODO offer option that if (multi)seasonal time series is provided
  # [msts()-object] then extract seasonalities
	wsize<- wsize #important as we modify S
  if(!is.null(extreg)){
    extreg <- as.matrix(extreg)
    if(is.null(colnames(extreg))){
      colnames(extreg) <- paste0("extreg_", 1:ncol(extreg))
    }
    # extreg containining missing is not recommended,
    # but algorithm will continue
    for(var.id in 1:dim(extreg)[2]){
      if(any(is.na(extreg[,var.id]))){
        extreg[,var.id]<- zoo::na.approx(extreg[,var.id], na.rm = FALSE)
        warning("extreg contains NAs - the decomposition may be inaccurate")
      }
    }
  }

  # sort S and get the indices of non-NA elements in x
  S <- sort(S)
	S<- S[S>1]
  idx <- which(!is.na(x))

  # derive the IC penalty factor
  if(ICpen=="HQC"){
    ICpen <- 2*log(log(length(idx)))
  } else if (ICpen=="AIC") {
    ICpen <- 2
  } else if (ICpen=="BIC"){
    ICpen <- log(length(idx))
  }

  # save the length of x
  n <- length(x)

  ## decide for use.autoregressive
  mmad <- function(z) stats::mad(z[z!=stats::median(z,na.rm=TRUE)],na.rm=TRUE)
  if(is.null(use.autoregressive)){
    if(mmad(diff(x)) < mmad(x)) use.autoregressive = TRUE else
      use.autoregressive = FALSE
  }

  #function for lagging the data
  get.lagged <- function(lag, Z){
    if(lag>=0) c(rep(NA, lag),  Z[(1+lag):(length(Z)) -lag ]) else
      c(Z[-(1:(-lag))], rep(NA, -lag) )
  }

  # function that calculates the running median
  get.rmed <- function(x, S, endrule="constant") {
    # x - numeric vector to be smoothed.
    # S - the seasonality.
    # endrule - character string indicating how the values at the beginning
    # and the end should be treated.

    # if seasonality is odd, then call directly the stats::runmed function
    # else calculate the average of running medians with S-1 and S+1
    if(as.logical(S%%2)){
      rmed <- as.numeric(stats::runmed(x, S, endrule=endrule))
    }  else{
      rmed <- 0.5*as.numeric(stats::runmed(x, S-1, endrule=endrule) +
                               stats::runmed(x, S+1, endrule=endrule))
    }  ## TODO end-rule could be replaced by something better, e.g. robets

    # replace the inner NAs by linear interpolation with zoo::na.approx over
    # the running median but leave the NAs on the bounds as the stats::runmed
    # does
    tmp <- rmed
    rmed[is.na(x)] <- NA
    rmed <- zoo::na.approx(rmed, na.rm = FALSE)
    tmpi <- which(is.na(rmed))
    rmed[tmpi] <- tmp[tmpi]

    return(rmed)
  }


  xw <- list() # list of weights
  xwl <- list() # list of weights including lags
  fit <- list() # list of fit
  fitl <- list() # list of fit on lags
  xreg <- list() # list of regressors
  xcode <- list() # list of coded regressor names

  # if use.trend == FALSE, then no running median is calculated and thus
  # no trend is estimated
  if(!use.trend){
    fit[[1]] <- numeric(n)
    xreg[[1]] <- cbind(rep(1,n), extreg)
    xcode[[1]] <- list(0) # 0 == intercept, -1== level,
    pnames <- c("intcpt", colnames(extreg), paste("S",S,sep=""))
  } else {
    rmed <- get.rmed(x, S=wsize)
    fit[[1]] <- rmed
    xreg[[1]] <- cbind(1,rmed , extreg)
    xcode[[1]] <- list(0,-1) # 0 == intercept, -1== level
    pnames <- c("intcpt", "trend", colnames(extreg), paste("S",S,sep=""))
  }

  if(!is.null(extreg)){
    for(i.col in seq_len(ncol(extreg))){
      xcode[[1]][[length(xcode[[1]])+1]] <- -1-i.col # 2,3,... == extreg
    }
  }



  ## function that calculates periodic-B-splines
  get.pbas.sparse <- function(n, period, dK = period / 4, ord = 4) {
    # n - length of splines to be calculated
    # period - seasonality
    # ord - order of the splines, e.g. 2 corresponds to quadratic, 3 to cubic
    # dK - equidistant distance
    # support will be 1:n

    stp <- dK
    x <- 1:period # must be sorted!
    lb <- x[1]
    ub <- x[length(x)]
    knots <- seq(lb - (0) * stp, ub + (1) * stp, by = stp)
    derivs <- numeric(length(x))

    ## some stuff adjusted from pbc-package to be faster
    degree <- ord - 1
    nKnots <- length(knots)
    Aknots <- c(knots[1] - knots[nKnots] + knots[nKnots - degree:1],
                knots, knots[nKnots] + knots[1:degree + 1] - knots[1])

    basisInterior <- splines::splineDesign(Aknots, x, ord, derivs, sparse=TRUE)
    basisInteriorLeft <- basisInterior[, 1:(ord - 1), drop = FALSE]
    basisInteriorRight <- basisInterior[, (ncol(basisInterior) - ord + 2
                                           ):ncol(basisInterior), drop = FALSE]
    basis <-cbind(basisInterior[, -c(
      1:(ord - 1), (ncol(basisInterior) - ord + 2):ncol(basisInterior)),
      drop = FALSE], basisInteriorLeft + basisInteriorRight)

    return(Matrix::Matrix(
      rep(Matrix::t(basis), length=n*as.numeric(dim(basis)[2])), n,
      dim(basis)[2], byrow=TRUE))
  }


  # function that calculates the functional periodic B-splines
  get.fpbas <- function(n, period=24, ord=4, period.is.integer=NULL, bexp = 2){
    # n - length of splines to be calculated
    # period - seasonality
    # ord - order of the splines, e.g. 2 corresponds to quadratic, 3 to cubic
    # bexp - basis of exponential grid
    # TODO at the moment only for integers, if not integer, increase Klist
    # but remove diag

    # estimate if smaller than 2 cycles then integer assumption is okay
    if(is.null(period.is.integer)){
      period.is.integer <- period %% 1 < period/n/2
    }
    if(period.is.integer){
      Klist <- bexp^seq(bexp, floor(log(period,bexp)-.5) )
    }  else Klist <- bexp^seq(log(ord, bexp), floor(log(period,bexp)) )

    kcum <- c(0, cumsum(Klist))
    i.k <- 1
    K <- Klist[i.k]

    BBl <- list()
    BBl[[1]] <- Matrix::sparseMatrix(i=seq_len(n),j=rep.int(1,n))+0

    for(i.k in seq_along(Klist)) {
			if( Klist[i.k] < ord ) next
      BBl[[1+i.k]] <- get.pbas.sparse(n, period =period,
                                      dK = period / Klist[i.k], ord = ord)
    }
    if(period.is.integer){
      BBl[[length(BBl)+1]] <-
        Matrix::sparseMatrix(i=seq_len(n),j=(seq_len(n)-1) %% period+1 ) +0
    }
    return(do.call("cbind", BBl))
  }

  # apply the get.fpbas function to each seasonality
  BBS <- list()
  for(i.S in seq_along(S)) {
#		if(S[i.S]>1)
		BBS[[i.S]] <- get.fpbas(n, S[i.S])
	}

  BBSk <- list()
  for(i.S in seq_along(S)) BBSk[[i.S]] <-
    textTinyR::sparse_Sums(BBS[[i.S]])/dim(BBS[[i.S]])[1]


  for(i.S in 1:max(length(S),1)){

    # vector of the remaining part, after subtraction of the trend
    xremainder <- x-fit[[i.S]]

    # calculate mad on the xremainder
    xmad <- mmad(xremainder)

    # calculate the weights for the estimation of weighted robust lasso
    # the higher the deviation from the median, the lower the weight
    xw[[i.S]] <- 1/( pmax(abs(xremainder/xmad) -K,0)+ 1)


    if(use.autoregressive){
      YMAT<-as.matrix(get.lagged(1,x))
      WMAT<-as.matrix(get.lagged(1,xw[[i.S]]))
      xwl[[i.S]]<- apply(cbind(WMAT, xw[[i.S]]), 1,min)
    } else {
      xwl[[i.S]]<- xw[[i.S]]
    }


    # prepare the list of all regressors for the decomposition
    BBSm <- list()

    Bki <- c()
    cnt <- 0
    xcode[[i.S+1]] <- list()
    xcode[[i.S+1]][[1]] <- 0

    # if trend is to be used, add the fitted running median
    lpen<- list()
    if(use.trend){
      BBSm[[length(BBSm)+1]] <- Matrix::Matrix(xreg[[i.S]][,2], sparse=TRUE)
      lpen[[length(lpen)+1]] <- 1
      xcode[[i.S+1]][[2]] <- -1
      Bki[length(Bki)+1] <- 1#dim(BBSm[[length(BBSm)]])[2]
    }

    if(!is.null(extreg)){
      for(i.col in seq_len(ncol(extreg))){
        BBSm[[length(BBSm)+1]] <- Matrix::Matrix(
          xreg[[i.S]][,1+use.trend+i.col], sparse=TRUE)
        lpen[[length(lpen)+1]] <- 1
        xcode[[i.S+1]][[length(xcode[[i.S+1]])+1]] <- -1-i.col
        Bki[length(Bki)+1] <- 1
      }
    }


    # prepare all remaining regressors and interactions for the estimation
		if(length(S)>0){
		  fi <- length(xcode[[i.S+1]] ) #1+use.trend# first include
		  for(i.SS in 1:i.S){
		    for(i.r in 1:dim(xreg[[i.S]])[2]){
		      if(all(xcode[[i.S]][[i.r]] < S[i.SS])){
		        xcode[[i.S+1]][[length(xcode[[i.S+1]])+1]] <-
		          c(xcode[[i.S]][[i.r]][xcode[[i.S]][[i.r]]!=0], S[i.SS])
		        if(fi < i.r){
		          BBSm[[length(BBSm)+1]] <- BBS[[i.SS]]*xreg[[i.S]][,i.r]
		          lpen[[length(lpen)+1]] <- 1/(BBSk[[i.SS]]+1/sqrt(n))
		          fi<- fi+1
		        } else {
		          BBSm[[length(BBSm)+1]] <- BBS[[i.SS]][,-1]*xreg[[i.S]][,i.r]
		          lpen[[length(lpen)+1]] <- 1/(BBSk[[i.SS]][-1]+1/sqrt(n))
		        }
		        Bki[length(Bki)+1] <- dim(BBSm[[length(BBSm)]])[2]
		      }
		    }
		  }
		}

    Bkicum <- c(0, cumsum(Bki))
    if(use.autoregressive) {
      BBSm[[length(BBSm)+1]] <- YMAT
      lpen[[length(lpen)+1]] <- 0
    }
    BB <- do.call("cbind", BBSm)

    idxa <- which(!is.na(xwl[[i.S]]))
    lpen <- unlist(lpen)
    lpen <- lpen/sum(lpen)

    # apply the weighted lasso
    model <- suppressWarnings(glmnet::glmnet(
      BB[idxa,],x[idxa], weights=xwl[[i.S]][idxa], penalty.factor=lpen,
      lambda=2^seq(-15,10,length=100), thresh=1e-7, dfmax=10*sqrt(length(idxa))
      ))


    # choose the best model IC-based
    RSS <- (1 - model$dev.ratio) * model$nulldev
    IC <- log(RSS) + model$df * ICpen / model$nobs*2
    idsel <- which.min(IC)
    xreg[[i.S+1]] <- array(,dim=c(dim(BB)[1], length(Bki)+1))

    # fill the xreg array with the estimated values
    xreg[[i.S+1]][,1] <- model$a0[idsel]
    for(i in 1:length(Bki)){
      itmpa <- (Bkicum[i]+1):Bkicum[i+1]
      if(length(itmpa)>1){
        xreg[[i.S+1]][,1+i] <-
          as.numeric( BB[,itmpa] %*% model$beta[itmpa, idsel] )
        # 1st in BBS is 1 (intercept)
      } else{
        xreg[[i.S+1]][,1+i] <-
          as.numeric( BB[,itmpa] * model$beta[itmpa, idsel] )
        # 1st in BBS is 1 (intercept)
      }
    }
    # fitted l-values
    if(use.autoregressive){
      fitl[[i.S+1]] <- apply(xreg[[i.S+1]],1,sum)
      attr(fitl[[i.S+1]], "phi") <- model$beta[dim(model$beta)[1], idsel]

      fit[[i.S+1]] <- fitl[[i.S+1]]*NA
      #fit[[i.S+1]][1] <- fitl[[i.S+1]][1]/(1- attr(fitl[[i.S+1]], "phi"))
#      fit[[i.S+1]][1] <- x[idx[1]]
			fit[[i.S+1]]<- stats::filter(c(x[idx[1]], fitl[[i.S+1]][-1]),
			                             model$beta[dim(model$beta)[1], idsel],
			                             method="recursive")
    } else {
      # fitted values
      fitl[[i.S+1]] <- numeric(0)
      fit[[i.S+1]] <- apply(xreg[[i.S+1]],1,sum)
    }
  }

  xremainder <- x - fit[[i.S+1]]
  xmad <- stats::mad(xremainder, na.rm=TRUE)

  # weights for next level approx, can be used e.g. in quantile regression
  # afterwards
  xw[[length(xw)+1]] <-  1/(pmax(abs(xremainder/xmad)-K, 0) + 1)

  # prepare the names of regressors
  pid <- (unique(unlist(xcode)))

  xregnames <- list()
  for(i.S in 1:(max(length(S),1)+1)){
    xregnames[[i.S]] <- numeric(0)
    for(i.n in seq_along(xcode[[i.S]])) xregnames[[i.S]][i.n] <-
        paste(pnames[match(xcode[[i.S]][[i.n]], pid)], collapse=":")
    dimnames(xreg[[i.S]]) <- list(NULL,xregnames[[i.S]])
  }

  results <- list(fit = fit[[length(fit)]], fitl = fitl[[length(fitl)]],
                  weights = xw[[length(xw)]],
                  components = xreg[[length(xreg)]])
  return(results)
}

#' Detects unreliable outliers in univariate time series data based on 
#' model-based clustering
#'
#' @description This function applies finite mixture modelling to compute
#' the probability of each observation being outliying data
#' in an univariate time series.
#' By utilizing the \link[mclust]{Mclust} package the data is
#' assigned in G clusters whereof one is modelled as an outlier cluster.
#' The clustering process is based on features, which are modelled to 
#' differentiate normal from outlying observation.Beside computing
#' the probability of each observation being outlying data also
#' the specific cause in terms of the responsible feature/ feature combination
#' can be provided.
#' @param data an one dimensional matrix or data frame without missing data;
#' each row is an observation.
#' @param S vector with numeric values for each seasonality present in data.
#' @param repetitions denotes the number of
#' repetitions to repeat the clustering.
#' By default set to 10. Allows to control the robustness and computational time
#' of the method.
#' @param proba denotes the threshold from which on an observation is considered
#' as being outlying data. By default is set to 0.5 (ranging from 0 to 1). Number of
#' outliers increases with decrease of proba threshold.
#' @param share controlls the size of the subsample used for estimation.
#' By default set to pmin(2*round(length(data)^(sqrt(2)/2)),
#' length(data))/length(data) (ranging from 0 to 1).
#' In combination with the repetitions parameter the
#' robustness and computational time of the method can be controlled.
#' @param decomp allows to perform seasonal decomposition on the original time series as pre-
#' processing step before feature modelling. By default set to TRUE.
#' @param PComp allows to use the principal components of the modelled feature matrix.
#' By default set to FALSE.
#' @param feat.inf logical value indicating whether influential features/ feature combinations
#' should be computed. By default set to FALSE.
#' @param out.par controls the number of artifially produced outliers to allow cluster
#' formation of oultier cluster. By default out.par ist set to 2. By increase it is assumed that
#' share of outliers in data increases. A priori it is assumed that
#' out.par * ceiling(sqrt(nrow(data.original))) number of observations are outlying observations.
#' @param detection.parameter denotes a parameter to regulate the
#' detection sensitivity. By default set to 1. It is assumed that the outlier cluster 
#' follows a (multivariate) Gaussian distribution parameterized by sample mean and a blown up
#' sample covariance matrix of the feature space. The covariance matrix is blown up
#' by detection.parameter * (2 * log(length(data)))^2.
#' By increase the more extrem outliers are detected.
#' @param max.cluster a single numeric value controlling the maximum
#' number of allowed clusters. By default set to 9.
#' @param G denotes the optimal number of clusters limited by the
#' max.cluster paramter. By default G is set to NULL and is automatically
#' calculated based on the BIC.
#' @param modelName denotes the geometric features of the covariance matrix.
#' i.e. "EII", "VII", "EEI", "EVI", "VEI", "VVI", etc.. By default modelName
#' is set to "VVV". The help file for \link[mclust]{mclustModelNames} describes
#' the available models. Choice of modelName influences the fit to the data as well as
#' the computational time.
#' @param ext.val denotes the number of observations for each side of an identified outlier,
#' which should also be treated as outliyng data. By default set to 1.
#' @param ... additional arguments for the \link[mclust]{Mclust} function.
#' @details The detection of outliers is addressed by
#' model based clustering based on parameterized finite Gaussian mixture models.
#' For cluster estimation the \link[mclust]{Mclust} function is applied.
#' Models are estimated by the EM algorithm initialized by hierarchical
#' model-based agglomerative clustering. The optimal model is selected
#' according to BIC.
#' The following features based on the introduced data are used in the clustering process:
#' \describe{
#'   \item{org.series}{denotes the scaled and potantially decomposed original time series.}
#'   \item{seasonality}{denotes determenistic seasonalities based on S.}
#'   \item{gradient}{denotes the summation of the two sided gradient of the org.series.}
#'   \item{abs.gradient}{denotes the summation of the absolute two sided gradient of
#'    org.series.}
#'   \item{rel.gradient}{denotes the summation of the two sided absolute gradient of the 
#'   org.series with sign based on left sided gradient in relation to the 
#'   rolling mean absolut deviation based on most relevant seasonality S.}
#'   \item{abs.seas.grad}{denotes the summation of the absolute two sided seasonal gradient of
#'   org.series based on seasonalties S.}
#'  }
#' In case PComp = TRUE, the features correspond to the principal components of the
#' introduced feature space.
#' @return a list containing the following elements:
#' \item{data}{numeric vector containing the original data.}
#' \item{outlier.pos}{a vector indicating the position of each outlier and the
#' corresponding neighboorhood controled by ext.val.}
#' \item{outlier.pos.raw}{a vector indicating the position of each outlier.}
#' \item{outlier.probs}{a vector containing all probabilities for each observation
#' being outlying data.}
#' \item{Repetitions}{provides a list for each repetition containing the estimated model,
#' the outlier cluster, the probabilities for each observation belonging to the estimated
#' clusters, the outlier position, the influence of each feature/ feature combination 
#' on the identified outyling data, and the corresponding probabilities 
#' after shift to the feature mean of each considered outlier, as well as the applied 
#' subset of the extended feature matrix for estimation (including artificially introduced 
#' outliers).
#' }
#' \item{features}{a matrix containg the feature matrix. Each column is a feature.}
#' \item{inf.feature.combinations}{a list containg the features/ feature comibinations,
#' which caused assignment to outlier cluster.}
#' \item{feature.inf.tab}{a matrix containing all possible feature combinations.}
#' \item{PC}{an object of class "princomp" containing the principal component analysis
#' of the feature matrix.}
#' @examples
#' \dontrun{
#' set.seed(1)
#' id <- 14000:17000
#' # Replace missing values
#' modelmd <- model_missing_data(data = GBload[id, -1], tau = 0.5,
#'                               S = c(48, 336), indices.to.fix = seq_len(nrow(GBload[id, ])),
#'                               consider.as.missing = 0, min.val = 0)
#' # Impute missing values
#' data.imputed <- impute_modelled_data(modelmd)
#' 
#' #Detect outliers
#' system.time(
#'   o.ident <- detect_outliers(data = data.imputed, S = c(48, 336))
#' )
#' 
#' # Plot of identified outliers in time series
#' outlier.vector <- rep(F,length(data.imputed))
#' outlier.vector[o.ident$outlier.pos] <- T
#' plot(data.imputed, type = "o", col=1 + 1 * outlier.vector,
#'      pch = 1 + 18 * outlier.vector)
#' 
#' # table of identified raw outliers and corresponding probs being outlying data
#' df <- data.frame(o.ident$outlier.pos.raw,unlist(o.ident$outlier.probs)[o.ident$outlier.pos.raw])
#' colnames(df) <- c("Outlier position", "Probability of being outlying data")
#' df
#' 
#' # Plot of feature matrix
#' plot.ts(o.ident$features, type = "o",
#'         col = 1 + outlier.vector,
#'         pch = 1 + 1 * outlier.vector)
#' 
#' # table of outliers and corresponding features/ feature combinations,
#' # which caused assignment to outlier cluster
#' # Detect outliers with feat.int = T
#' set.seed(1)
#' system.time(
#'   o.ident <- detect_outliers(data = data.imputed, S = c(48, 336), feat.inf = T)
#' )
#' feature.imp <- unlist(lapply(o.ident$inf.feature.combinations,
#'                              function(x) paste(o.ident$feature.inf.tab[x], collapse = " | ")))
#' 
#' df <- data.frame(o.ident$outlier.pos.raw,o.ident$outlier.probs[o.ident$outlier.pos.raw],
#'                  feature.imp[as.numeric(names(feature.imp)) %in% o.ident$outlier.pos.raw])
#' colnames(df) <- c("Outlier position", "Probability being outlying data", "Responsible features")
#' View(df)
#' }
#' @export
#' @seealso \code{\link[tsrobprep]{model_missing_data}},
#' \link[tsrobprep]{impute_modelled_data},
#' \link[tsrobprep]{auto_data_cleaning}
#' @importFrom mclust mclustBIC
#'
# detect outliers function - modified
detect_outliers <- function(data,
                            S,
                            proba = 0.5,
                            share = NULL,
                            repetitions = 10,
                            decomp = T,
                            PComp = F,
                            detection.parameter = 1,
                            out.par = 2,
                            max.cluster = 9,
                            G = NULL,
                            modelName = "VVV",
                            feat.inf = F,
                            ext.val = 1,
                            ...) {
  # transform data into matrix format
  data.original <- as.matrix(data)

  # validate the variables' correctness - basic validation
  if(sum(is.na(data.original)) != 0) stop("Please impute missing data before usage.")
  if(is.null(share)){
    subset.length <- pmin(2*round(length(data.original)^(sqrt(2)/2)),
                                  length(data.original))
    # Set share to 1 if number of observations are smaller than 1000
    if(subset.length < 1000){
      share <- 1
    } else {
      share <- subset.length/length(data.original)
    }
  }
  
  if(is.null(G)){ # Create grid for search of clusters
    phi<- (sqrt(5)+1)/2
    G<- round((phi^(2+1:log(max.cluster-1, phi)))/sqrt(5) )
  }
  
  ## TODO Disable seasonal features if data is too short
  if(sum(nrow(data.original) <= S) < length(S))  S <- S[which(S<nrow(data.original))]

  ## 1.) time series decomposition
  if (decomp == T){
    # Transform data to ts/msts object
    #msts_data <- msts(data.original, seasonal.periods = S, ts.frequency = S[1])
    ts_data <- stats::ts(data.original[,1], frequency = S[1])
    # Seasonal decompostion with stl()
    # TODO robust = TRUE allows robust fitting
    dts_data <- stats::stl(ts_data, s.window = "periodic", robust = TRUE)
    # plot(dts_data$time.series)

    # Use remainder of stl object
    data.original <- as.matrix(dts_data$time.series[,3])
  }

  ## 2.) Feature engineering
  # Original sereis
  org.series <- as.numeric(data.original)
  # One sided gradients
  dyl <- c(0, diff(org.series)) # Differences (y_{i+1} - y_{y-i})
  dyr <- c(diff(org.series), 0) # Differences (y_{i+1} - y_{y-i})
  # Summation of left and right sided gradient
  gradient <- dyl + dyr # Summation of dyl and dyr
  ## Absolute two sided gradient
  abs.gradient <- abs(dyl) + abs(dyr) - abs(dyl + dyr) # Absolute gradient
  # Summation of dyl and dyr in accordance to sign(dyl)
  pos.neg.gradient <- dyl + sign(dyl) * abs(dyr) 
  # Computation of seasonal mad matrix based on strongest seasonality S 
  seas.mat <- matrix(, S[1], ceiling(length(pos.neg.gradient)/S[1]))
  seas.mat[1:length(pos.neg.gradient)]<- pos.neg.gradient
  seas.vol <- rep_len(apply(seas.mat, 1, stats::mad, na.rm=TRUE),
                      length.out=length(pos.neg.gradient))
  # TODO Tuning parameter to control smoothness of rolling median of seasonal mad
  H<-1
  seas.vol.ma <- sapply(1:length(seas.vol), function(x)
    stats::median(seas.vol[pmax(x-H,1):pmin(x+H,length(seas.vol))]))
  # Add noise to aviod zero pos.neg.gradient
  pos.neg.gradient_sc <- scale(pos.neg.gradient + 10^-100) 
  # Relative gradient
  # Add noise to aviod zero pos.neg.gradient
  rel.gradient <- as.numeric(pos.neg.gradient_sc / (seas.vol.ma + 10^-100))
  # Seasonal features
  names <- vector()
  seas_feat_mat <- matrix(0, nrow = nrow(data.original), ncol = 0)
  for (i in seq_len(length(S))) {
    # Linear trend
    seas_feat_mat <- cbind(seas_feat_mat,
                           as.numeric(scale(seq_along(org.series) %% S[i])))
    # Seasonal differences (y_{i+S} - y_{y-i})
    dylS <- c(rep(0, S[i]), diff(org.series, S[i]))
    # Seasonal differences (y_{i+S} - y_{y-i})
    dyrS <- c(diff(org.series, S[i]), rep(0, S[i]))
    #Absolute seasonal gradient
    seas_feat_mat <- cbind(seas_feat_mat,
                           abs(dylS) + abs(dyrS) - abs(dylS + dyrS))
    # Colnames
    names <- c(names, paste0("seasonality.", S[i]), paste0("abs.seas.grad.S", S[i]))
  }

  # Merge features into matrix
  feat.internal <- cbind(org.series, gradient, abs.gradient, rel.gradient, seas_feat_mat)
  colnames(feat.internal)[5:ncol(feat.internal)] <- names
  
  # TODO Allow users to choose features or to introduce new features

  # scale features
  feat.internal.scaled <- apply(feat.internal,2,scale)

  ## 3.) Principal component analysis
  if(PComp == T){
    system.time(
      pf<- stats::princomp(feat.internal.scaled)
    )
    # sdev = standard deviations of PCs
    # scores of the supplied data on the PCs
    feature.mat<- t(pf$sdev * t(pf$scores))
    # Restore principal components
    PCs <- pf

  } else {
    feature.mat <- feat.internal.scaled # scaled but original features
    PCs <- NULL
  }

  ## 4.) Expand outlier cluster by artificially introducing
  # outliers based on blown up covariance matrix
  mu <- apply(feature.mat, 2, mean) # Calculation of feature means
  sigma <- stats::cov(feature.mat) # Calculation of covariance matrix

  # Draw artificial outliers from outlier distribution
  size_out_cluster <- out.par * ceiling(sqrt(nrow(data.original))) # to extend outlier cluster size relativ for small cluster

  # Creation of outlier cluster
  pts <- MASS::mvrnorm(size_out_cluster, mu,
                       sigma * detection.parameter * (2 * log(nrow(data.original)))^2)
  # add small noise
  # TODO 10^-5 denotes noise parameter. May become tuning parameter
  # in case user can introduce features.
  zz <- feature.mat + stats::rnorm(length(feature.mat)) * 10^-5

  ## Repetitions
  REP.list <- list()
  time.feature.inf.list <- list()
  feature.list <- list()
  for(j in 1:repetitions){
    # 5.) Subsetting of data
    zz_ss_ind <- sample(seq_len(nrow(data.original)),
                        ceiling(share * nrow(data.original)), replace = FALSE)
    pts_ss_ind <- sample(seq_len(nrow(pts)),
                         out.par * ceiling(sqrt(share * nrow(data.original))), replace = F)
    zzaug <- rbind(zz[zz_ss_ind,], pts[pts_ss_ind,])

    # Partitioning: Assigning artificial outlier cluster into one group
    pp<- c(seq_len(length(zz_ss_ind)), rep(length(zz_ss_ind)+1, length(pts_ss_ind)))

    ## 6.) Run Mclust algorithm
    mod <- mclust::Mclust(zzaug, G = G, modelName = modelName, verbose = FALSE, partition = pp)

    # Extract outlier cluster
    # Obtain variances
    vars <- summary(mod, parameters = TRUE)$var
    # Function to extract diagonal of covariance matrix
    get.dvar <- function(z) diag(vars[, , z])
    # Extract diagonal of covariance matrix for each cluster
    dvars <- sapply(1:dim(vars)[3], get.dvar)
    # Receive standard deviation of each cluster
    dec <- apply(sqrt(dvars), 2, sum)

    dens.mat <- matrix(, nrow = nrow(data.original), ncol = length(dec))
    # Compute densities of each observation to belong to clusters
    for (g in seq_len(length(dec))) {
      # Mean of cluster g
      mu <- mod$parameters$mean[, g]
      # Covariance of cluster g
      sigma <- mod$parameters$variance$sigma[, , g]
      # Compute densities of each observation to belong to clusters
      dens.mat[, g] <- mclust::dmvnorm(zz, mu, sigma)
    }

    ## 7.) Compute probability of outliers belonging to kth cluster
    probs <- t(apply(dens.mat, 1, function(x) x / sum(x)))

    # Logical values inidcating whether value is an outlier or not
    outlier.vec<-rep(F,nrow(data.original))
    outlier.vec[which(probs[,which.max(dec)] >= proba)] <- TRUE

    ## 8.) Compute feature influence
    # compute feauture combinations
    if(sum(outlier.vec) != 0 & feat.inf == T){
      # Grid of all combinations based on all features
      exp.gri_base <- expand.grid(replicate(ncol(feature.mat), c(0,1), simplify=FALSE))[-1,]
      # Set deterministic features equal to zero

      exp.gri_base[,seq(5,ncol(feature.mat),2)[1:length(S)]] <- 0
      # Use only those feautre combinations where mean shift 

      # is reasonable...
      exp.gri <- unique(exp.gri_base)
      rownames(exp.gri) <- 1:nrow(exp.gri)
      colnames(exp.gri) <- sapply(seq_len(ncol(exp.gri)), function(x) colnames(feature.mat)[x])
      feat.inf.mat <- matrix(F, nrow = nrow(data.original), ncol = nrow(exp.gri)) # matrix to restore influence of features on outlier assignment
      feat.inf.grid.mat <- feat.inf.mat[outlier.vec,]
      probs.outliers.mean.s <- array(NA, dim = c(nrow(data.original), nrow(exp.gri), mod$G)) # Array to restore probabilities of observations to belong to outlier cluster after feautre manipulation
      outlier.feat.list <- list() # list of influential feature combination
      for(k in 1:nrow(exp.gri)){
        inde <- which(exp.gri[k,]==1)
        # Shift feature values to feature means of identified outliers
        feature.mat.adj <- feature.mat
        # Shift mean of all outliers and for the corresponding feature combination at once
        # TODO How to calculate mean?
        feature.mat.adj[outlier.vec, inde] <-  matrix(rep(apply(as.matrix(feature.mat[-outlier.vec, inde]),2, mean), sum(outlier.vec)),
                                                      nrow = sum(outlier.vec), ncol = length(inde), byrow =T)

        # Objects to restore densities of manipulated outliers
        dens.mat.mani <- matrix(,nrow=sum(outlier.vec), ncol=length(dec))

        # Compute densities for manipulated outliers belonging to clusters
        for(g in 1:length(dec)){
          # Receive means and sigma for specific cluster
          mu <- mod$parameters$mean[,g]
          sigma <- mod$parameters$variance$sigma[,,g]

          # Compute densities for manipulated outlier vector
          dens.mat.mani[,g]<-mclust::dmvnorm(feature.mat.adj[outlier.vec, ,drop=FALSE], mu, sigma)
        }

        # Compute probability of outliers to belong to kth cluster
        probs.outliers.mean.s[which(outlier.vec==T),k,] <- t(apply(dens.mat.mani, 1, function(x) x/sum(x)))

        # Assign logical value to specific feature combination which caused outlier assignment
        moved.outliers <- which(probs.outliers.mean.s[which(outlier.vec==T),
                                                      k, which.max(dec)] <= proba)
        feat.inf.mat[which(outlier.vec==T)[moved.outliers], k] <- T
      }

      ## 9) Create output list (reduced version)
      # nrow = identified outliers of repetition i
      # ncol = feature combinations
      feat.inf.grid.mat <- feat.inf.mat[outlier.vec, ,drop=FALSE]

      # Select only features which caused class shift
      outlier.feat.list <- sapply(as.character(which(outlier.vec == T)),function(x) NULL)
      for (jj in 1:sum(outlier.vec)){
        # matrix containg only feature combination which caused class shift
        exp.gri.mat <- exp.gri[which(feat.inf.grid.mat[jj, ,drop=FALSE] == T),]

        ## Compute all combinations
        r.sums <- rowSums(exp.gri.mat) # how often occurs each feature combination

        l.el <- list()
        for (l in seq_len(ncol(feature.mat))){
          l.el <- c(l.el,sapply(which(r.sums == l), function(x) which(exp.gri.mat[x,] == 1), simplify = F))
        }

        name <- paste(which(outlier.vec==T)[jj])
        outlier.feat.list[[name]] <- l.el

        ind <- vector()
        for (jjj in 1:length(l.el)){
          ind <- c(ind, which(lapply(l.el, function(x) sum(l.el[[jjj]] %in% x) == length(l.el[[jjj]])) ==T))
        }

        # Feature combination, which caused outlier assignment - unqiue
        outlier.feat.list[[name]] <- as.numeric(names(l.el[which(table(ind)==1)]))
      }
    } else {
      probs.outliers.mean.s <- NULL
      outlier.feat.list <- NULL

    }

    ## 10) Restore subsample results
    REP.list[[paste("Rep",j)]] <- list("Model" = mod,
                                       "Outlier.clus" = which.max(dec),
                                       "Probs" = probs,
                                       "Outlier.pos" = outlier.vec,
                                       "outlier.feat.list" = outlier.feat.list,
                                       "ProbsMeanShift" = probs.outliers.mean.s,
                                       "Subset" = zzaug
    )
  }

  ## 11) Average results
  # All detected outliers over all repetitions based on proba
  outlier.vec.rep <- unique(unlist(lapply(REP.list, function(x) which(x$Outlier.pos==T))))
  feature.list <- sapply(as.character(outlier.vec.rep),function(x) NULL)

  # Restore Influential Features over all repetitions and for each detected outlier based on proba
  if( length(outlier.vec.rep ) != 0 & feat.inf == T){
    for( k in names(feature.list)){
      for( kk in seq_len(repetitions)){
        # unique = drop duplicated feature combinations over repetitions
        feature.list[[k]] <- unique(c(feature.list[[k]],
                                        unlist(REP.list[[kk]]$outlier.feat.list[which(names(REP.list[[kk]]$outlier.feat.list) == k)])))
      }
    }
  }

  # Probabilities of being outlying data for each observation and over all repetitions
  probs.mat <- matrix(unlist((lapply(REP.list, function(x) x$Probs[,x$Outlier.clus]))), ncol = length(REP.list), byrow = F)

  ## Results
  # outlier.pos.raw
  probs.outlier.av <- apply(probs.mat,1,mean)
  outlier.pos.av <- which(probs.outlier.av > proba) # Final outliers based on averaging and above proba threshold

  #outlier.pos augmented 
  outlier.pos.aug <- as.vector(unique(sapply(outlier.pos.av, function(x) pmax(x-ext.val,1):pmin(x+ext.val, nrow(data.original)))))
  if(length(outlier.pos.aug)==0) outlier.pos.aug <- NULL
  
  # feature grid
  if( length(outlier.vec.rep ) != 0 & feat.inf == T){
    exp.grid <- apply(exp.gri,1,function(x) names(which(x == 1)))
  } else {
    exp.grid <- NULL
  }

  result <- list(data = data,
                 outlier.pos = outlier.pos.aug,
                 outlier.pos.raw = outlier.pos.av,
                 outlier.probs = probs.outlier.av,
                 repetitions =  REP.list,
                 features = feature.mat,
                 inf.feature.combinations = feature.list,
                 feature.inf.tab = exp.grid,
                 pc = PCs
  )
  return(result)
}


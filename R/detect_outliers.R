#' Detects unreliable outliers in univariate time series data based on finite
#' mixture modelling
#'
#' @description This function applies finite mixture modelling to compute
#' the probability of each observation being an outlier
#' in an univariate time series.
#' By utilizing the \link[mclust]{Mclust} package the data is
#' assigned in k clusters whereof one contains outliers.
#' The clustering process is based on features, which are
#' specifically desinged for outlier detection in time series data.
#' @param data an one dimensional matrix or data frame without missing data;
#' each row is an observation.
#' @param S vector with numeric values for each seasonality present in data.
#' @param feat.int a list of logical values indicating which features should
#' be applied in clustering algorithm.
#' \describe{
#'   \item{org.s}{denotes the scaled original time series.}
#'   \item{lin.tr}{denotes linear trends based on seasonalities S.}
#'   \item{grad}{denotes the gradient of scaled original time series.}
#'   \item{abs.grad}{denotes absolute gradient of
#'    scaled original time series.}
#'   \item{abs.seas.grad}{denotes the absolute seasonal gradient of
#'   scaled original time series based on seasonalties S.}
#'}
#' @param model.select.iter denotes the number of
#' iterations to find the optimal number of clusters as well as
#' the optimal model of the covariance matrix.
#' By default set to 10. Results become more robust
#' with increasing iterations but lead
#' likewise to increasing computational time.
#' @param outlier.detect.iter denotes the number of iterations
#' for outlier detection based on variing subsamples
#' controlled by the share parameter. By default set to 10.
#' Results become more robust with increasing iterations but lead
#' likewise to increasing computational time.
#' @param proba denotes threshold ranging from 0 to 1 from which an
#' observation is considered being an outlier. Number of outlier increases
#' with decreasing threshold. By default set to 0.1, implying that in average
#' an observation obtains a probability in each iteration of 0.1 belonging to
#' the outlier cluster.
#' @param share controlls the size of the subsample used for estimation.
#' By default set to 0.3 (ranging from 0 to 1). Controlls the computational time
#' and the robustness of the method.
#' @param detection.parameter denotes a parameter to regulate the
#' detection sensisitivity. By default set to 1. The smaller the more outliers.
#' @param out.par controls the number of artifially produced outliers in
#' relation to the subsample size controlled by the share parameter. By default
#' out.par ist set to 2. By increase a priori it is assumed that
#' share of outliers in data increases.
#' @param noise.par controls strenght of noise added to feature matrix to avoid
#' zero variance issues by applying bivariate features. By default set to 10^-5,
#' strength of noise decreases with decreasing noise.par parameter.
#' @param mar denotes a margin controlling the number of adjacent
#' values around an identified outlier which are likewise considered
#' as outliers. By default set to 1 so that the one most closest neighbours
#' of an identified outlier on each side are also treated as outliers.
#' @param max.cluster a single numeric value controlling the maximum
#' number of clusters allowed. By default set to 9.
#' @param G denotes the optimal number of clusters limited by the
#' max.cluster paramter. By default G is set to NULL and automatically
#' calculated based on the BIC.
#' @param modelName denotes the geometric features of the covariance matrix.
#' i.e. "EII", "VII", "EEI", "EVI", "VEI", "VVI", etc.. By default modelName
#' is set to NUll and automatically calcualted based on BIC.
#' The help file for \link[mclust]{mclustModelNames} describes
#' the available models.
#' @param ... additional arguments for the \link[mclust]{Mclust} function.
#' @details The detection of outliers is addressed by
#' model based clustering based on parameterized finite Gaussian mixture models.
#' For cluster estimation the \link[mclust]{Mclust} function is applied.
#' Models are estimated by the EM algorithm initialized by hierarchical
#' model-based agglomerative clustering. The optimal model can be selected
#' according to BIC.
#' @return An object of class "tsrobprep" containing the following elements:
#' \describe{
#' \item{original.data}{an numeric vector containing the original data.}
#' \item{outlier.probs}{an numeric vector containing the averaged probability}
#' \item{outlier.probs.mat}{a matrix containg the probability for each iteration
#' that observation is belonging to the outlier cluster.
#' Each row is an observation and each column an iteration.}
#' \item{outlier.pos}{a logical vector indicating the position of each outlier.}
#' \item{outlier.pos.aug}{a logical vector indicating the position of
#' each outlier including neigbouring values based on the mar parameter.}
#' \item{estimated.models}{a list containing each estimated model.}
#' \item{BIC}{an mclustBIC object containing the
#' Bayesian Information Criterion for the specified mixture models numbers of
#' clusters. Auxiliary information returned as attributes.}
#' }
#' @examples
#' \dontrun{
#' set.seed(1)
#' id <- 14000:17000
#' # Replace missing values
#' modelmd <- model_missing_data(data = GBload[id, -1], tau = 0.5,
#'  S = c(48, 336), indices.to.fix = seq_len(nrow(GBload[id, ])),
#'  consider.as.missing = 0, min.val = 0)
#' # Impute missing values
#'  data.imputed <- impute_modelled_data(modelmd)
#' # Detect outliers
#' o.ident <- detect_outliers(data = data.imputed, S = c(48, 336),
#'                            model.select.iter = 1,
#'                            outlier.detect.iter = 1)
#' # Plot of identified outliers in time series
#' plot(data.imputed, type = "o", col=1 + o.ident$outlier.pos.aug,
#'      pch = 1 + 18 * o.ident$outlier.probs)
#'
#' # Plot of feature matrix
#' plot.ts(o.ident$features, type = "o",
#'         col = 1 + o.ident$outlier.pos,
#'         pch = 1 + 18 * o.ident$outlier.probs)
#' }
#' @export
#' @seealso \link[tsrobprep]{model_missing_data},
#' \link[tsrobprep]{impute_modelled_data},
#' \link[tsrobprep]{auto_data_cleaning}
#' @importFrom mclust mclustBIC
#'
detect_outliers <- function(data,
                            S,
                            model.select.iter = 10,
                            outlier.detect.iter = 10,
                            proba = 0.1,
                            share = 0.3,
                            detection.parameter = 1,
                            out.par = 2,
                            noise.par = 10^-5,
                            mar = 1,
                            max.cluster = 9,
                            G = NULL,
                            modelName = NULL,
                            feat.int = list(org.s = TRUE,
                                            grad = TRUE,
                                            abs.grad = TRUE,
                                            abs.seas.grad  = TRUE,
                                            lin.tr = TRUE),
                            ...) {
  data.original <- as.matrix(data) # matrix format

  ###validate the variables' correctness - basic validation
  if (dim(data.original)[1] == 0) {
    stop("Provided data is of length 0.")
  } else nobs <- dim(data.original)[1]
  if (length(S) == 0 | any(sapply(S, function(x) x %% 1 == 0) == FALSE))
    stop("No or incorrect seasonalities provided")
  if (!is.numeric(model.select.iter)) stop("Provided runs are not numeric.")
  if (!is.numeric(outlier.detect.iter)) stop("Provided runs are not numeric.")
  if (!is.numeric(proba) | proba > 1 | proba < 0)
    stop("Provided probability is incorrect.")
  if (!is.numeric(detection.parameter))
    stop("Provided detection.parameter is incorrect.")
  if (!is.null(G) & !is.numeric(G)) stop("Provide numeric value.")
  if (!is.null(modelName) & !is.character(modelName))
    stop("Provided modelName is incorrect.")
  if (!is.numeric(out.par)) stop("Provided out.par is not numeric.")
  if (!is.numeric(noise.par)) stop("Provided out.par is not numeric.")
  if (!is.null(G) & !is.null(modelName)) model.select.iter <- 0
  if (sum(is.na(data.original) > 0)) stop("Remove missing values.")
  if ((! TRUE %in% feat.int)) stop("No features provided.")
  ## Feature engineering
  # Original sereis
  org.series <- as.numeric(scale(data.original)) # scale data
  # Gradient
  dyl <- c(0, diff(org.series)) # Differences (y_{i+1} - y_{y-i})
  dyr <- c(diff(org.series), 0) # Differences (y_{i+1} - y_{y-i})
  gradient <- dyl + dyr # Summation of dyl and dyr
  # Absolute gradient
  abs.gradient <- abs(dyl) + abs(dyr) - abs(dyl + dyr) # Absolute gradient
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
    names <- c(names, paste0("lin.trend.S", S[i]), paste0("abs.seas.grad.S", S[i]))
  }

  # Merge features into matrix
  feat.internal <- cbind(org.series, gradient, abs.gradient, seas_feat_mat)
  colnames(feat.internal)[4:ncol(feat.internal)] <- names

  # Select relevant features based on user input
  rel_feat_ind <- c(unlist(feat.int)[1:3],
                    rep(unlist(feat.int)[4:5], length(S)))
  feature.mat <- feat.internal[, rel_feat_ind]

  # Check whether proved number of observations is sufficent large
  if (max.cluster * ncol(feature.mat) > share * nobs)
    stop("provided number of observations are too small.")

  ## Expand outlier cluster
  mu <- apply(feature.mat, 2, mean) # Calculation of feature means
  sigma <- stats::cov(feature.mat) # Calculation of covariance matrix

  # Draw artificial outliers from outlier distribution
  size_out_cluster <- out.par * ceiling(sqrt(nobs))

  # Production of nn samples from specified
  # multivariate normal distribution (mu, sigma, K)
  pts <- MASS::mvrnorm(size_out_cluster, mu,
                       sigma * detection.parameter * (2 * log(length(data)))^2)
  # add small noise
  zz <- feature.mat + stats::rnorm(length(feature.mat)) * noise.par

  # Objects to restore results
  models <- list() # Restore estimated models
  # Determine optimal number of clusters and
  # model structure of covariance matrix
  if (is.null(G) |  is.null(modelName)) {
    BIC <- NULL
    for (i in 1:model.select.iter) {
      # Subsample to increase run time
      zz_ss_ind <- sample(seq_len(nobs),
                          ceiling(share * nobs), replace = FALSE)
      out_aug_ind <- sample(seq_len(nrow(pts)),
                            out.par * ceiling(sqrt(share * nobs)), replace = F)
      zzaug <- rbind(zz[zz_ss_ind, ], pts[out_aug_ind, ])
      models[[i]] <- mclust::Mclust(zzaug,
                                    G = 1 : max.cluster, verbose = FALSE)
      BIC <- mclust::mclustBICupdate(BIC, models[[i]]$BIC)
    }
    # Assign optimal number of clusters and model structure
    # of covariance matrix to variables
    G <- which.max(apply(BIC, 1, function(x) max(stats::na.omit(x))))
    modelName <- colnames(BIC)[which.max(apply(BIC, 2,
                                               function(x) max(stats::na.omit(x))))]
  }
  # Objects to restore results
  outlier.vec <- rep(F, nobs) # Vector to restore positions of outliers
  # Matrix to restore probabilities that observations belongs to outlier cluster
  probs.out.cluster <- matrix(NA, nrow = nobs,
                              ncol = outlier.detect.iter + model.select.iter)

  for (j in seq_len(outlier.detect.iter + model.select.iter)) {
    if (j > model.select.iter) {
      while (length(models) < j) {
        # Subsample to increase run time
        zz_ss_ind <- sample(seq_len(nobs), ceiling(share * nobs), replace = F)
        out_aug_ind <- sample(seq_len(nrow(pts)),
                              out.par * ceiling(sqrt(share * nobs)),
                              replace = F)
        zzaug <- rbind(zz[zz_ss_ind, ], pts[out_aug_ind, ])
        ### TODO ###
        # Make out why mclust algortihm crahses sometimes --> not clear yet
        # Solution by now: Reestimation

        # G = optimal number of mixture components
        # modelName = model to be fit in EM pahse
        models[[j]] <- mclust::Mclust(zzaug, G = G, modelNames = modelName, verbose = FALSE)
        if (length(models) < j)  warning("Mclust provided empty model,
        which led to reestimation with new subset of data")
      }
    }
    # Extracting cluster with highest variance (Outlier cluster)
    vars <- summary(models[[j]], parameters = TRUE)$var
    # Function to extract diagonal of covariance matrix
    get.dvar <- function(z) diag(vars[, , z])
    # Extract diagonal of covariance matrix for each cluster
    dvars <- sapply(1:dim(vars)[3], get.dvar)
    # Receive standard deviation of each cluster
    dec <- apply(sqrt(dvars), 2, sum)
    # Matrix to store density of each observation to belong to clusters
    dens.mat <- matrix(, nrow = nobs, ncol = length(dec))
    # Compute densities of each observation to belong to clusters
    for (g in seq_len(length(dec))) {
      # Mean of cluster g
      mu <- models[[j]]$parameters$mean[, g]
      # Covariance of cluster g
      sigma <- models[[j]]$parameters$variance$sigma[, , g]
      # Compute densities of each observation to belong to clusters
      dens.mat[, g] <- mclust::dmvnorm(zz, mu, sigma)
    }
    # Compute probability of outliers to belong to kth cluster
    probs <- t(apply(dens.mat, 1, function(x) x / sum(x)))
    # Restore only probabilities of outlier cluster
    probs.out.cluster[, j] <- probs[, which.max(dec)]
  }
  # Vector with numeric values indicating
  # whether observation is an outlier or not
  probs_out_vec <- apply(probs.out.cluster, 1, mean)
  # Logical values inidcating whether value is an outlier or not
  outlier.vec[which(probs_out_vec >= proba)] <- TRUE
  # Vector with logical values indicating whether
  # observation (with adjacent values) is an outlier or not
  outlier.aug <- rep(FALSE, nobs)
  outlier.aug[unique(unlist(sapply(which(outlier.vec == TRUE),
                                   function(x) pmax(x - mar, 1) : pmin(x + mar, nobs))))] <- TRUE
  result <- list(original.data = data.original,
                 outlier.probs =  probs_out_vec,
                 estimated.models = models,
                 BIC = BIC,
                 features = feature.mat,
                 outlier.pos = outlier.vec,
                 outlier.pos.aug = outlier.aug,
                 outlier.probs.mat = probs.out.cluster)
  #class(result) <- "tsrobprep"
  return(result)
}

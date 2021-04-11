#' Impute modelled missing time series data
#'
#' Returns a matrix or a list of matrices with imputed missing values or
#' outliers. As argument the function requires an object of class "tsrobprep"
#' and the quantiles to be imputed.
#' @param object an object of class "tsrobprep" that is an output of function
#' model_missing_data.
#' @param tau the quantile(s) of the missing values to be imputed. tau should be
#' a subset of the quantile values present in the "tsrobprep" object. By default
#' all quantiles present in the object are used.
#' @return A matrix or a list of matrices with imputed missing values or
#' outliers.
#' @examples
#' model.miss <- model_missing_data(data = GBload[,-1], S = c(48,7*48), tau = 0.5,
#'               no.of.last.indices.to.fix = dim(GBload)[1], consider.as.missing = 0,
#'               min.val = 0)
#' model.miss$estimated.models
#' model.miss$replaced.indices
#' new.GBload <- impute_modelled_data(model.miss)
#' @export
#' @seealso \link[tsrobprep]{model_missing_data}, \link[tsrobprep]{detect_outliers},
#' \link[tsrobprep]{auto_data_cleaning}
impute_modelled_data <- function(object, tau = NULL){

  if(class(object)!="tsrobprep") stop("provided object is not of tsrobprep class")
  if(!all(tau %in% object$tau)) stop("provided tau is incorrect")

  if(length(object$replaced.values) == 0){
    warning("nothing to impute")
    return(object$original.data)
  }

  if(is.null(tau)) tau <- object$tau


  matrix.list <- list()

  for(tau.no in seq_along(tau)){

    object.tau.no <- match(tau[tau.no], object$tau)

    tau.char <- as.character(tau[tau.no])
    matrix.list[[tau.char]] <- object$original.data

    for(var.id in seq_len(dim(object$original.data)[2])){
      missing.index <- unlist(object$replaced.indices[[var.id]])
      matrix.list[[tau.char]][missing.index, var.id] <- object$replaced.values[[var.id]][missing.index,tau.no]
    }
  }
  if(length(matrix.list) == 1){
    return(matrix.list[[1]])
  } else return(matrix.list)
}

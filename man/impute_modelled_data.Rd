% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute_modelled_data.R
\name{impute_modelled_data}
\alias{impute_modelled_data}
\title{Impute modelled missing time series data}
\usage{
impute_modelled_data(object, tau = NULL)
}
\arguments{
\item{object}{an object of class "tsrobprep" that is an output of function
model_missing_data.}

\item{tau}{the quantile(s) of the missing values to be imputed. tau should
be a subset of the quantile values present in the "tsrobprep" object. By
default all quantiles present in the object are used.}
}
\value{
A matrix or a list of matrices with imputed missing values or
outliers. \insertNoCite{*}{tsrobprep}
}
\description{
Returns a matrix or a list of matrices with imputed missing values or
outliers. As argument the function requires an object of class "tsrobprep"
and the quantiles to be imputed.
}
\examples{
\dontrun{
model.miss <- model_missing_data(
    data = GBload[,-1], S = c(48,7*48),
    no.of.last.indices.to.fix = dim(GBload)[1], consider.as.missing = 0,
    min.val = 0
)
model.miss$estimated.models
model.miss$replaced.indices
new.GBload <- impute_modelled_data(model.miss)
}
}
\references{
\insertAllCited{}
}
\seealso{
\link[tsrobprep]{model_missing_data},
\link[tsrobprep]{detect_outliers}, \link[tsrobprep]{auto_data_cleaning}
}

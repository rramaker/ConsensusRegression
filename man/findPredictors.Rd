% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/findPredictors.R
\name{findPredictors}
\alias{findPredictors}
\title{An implementation of the step function to find covariates that might be important to include in regression}
\usage{
findPredictors(x, response, method, direction = "both")
}
\arguments{
\item{x}{List generated by the getData function}

\item{response}{Column header your response of interest in info file (i.e."Disease")}

\item{method}{Intended regression method.  Supports all base glm methods. (i.e. "gaussian", "binomial", "poisson" etc.)}

\item{direction}{Direction for step function.  Can be "backwards", "forwards", or "both".  Defaults to both}
}
\value{
Returns a vector of infoFile variables likely to be important in a model
}
\description{
An implementation of the step function to find covariates that might be important to include in regression
}


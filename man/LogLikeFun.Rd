% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunsForOptimalV2.R
\name{LogLikeFun}
\alias{LogLikeFun}
\title{LogLikeFun}
\usage{
LogLikeFun(bb, ixx, iyy, iw, iZZ)
}
\arguments{
\item{bb}{initial values for the intercept and slope coefficients}

\item{ixx}{continuous predictor}

\item{iyy}{binary outcome}

\item{iw}{the weighted parameter}

\item{iZZ}{covariates to be incorporated in the model}
}
\value{
a single log-likelihood value
}
\description{
caculate the log-likelihood value for Maximum Likelihood Estimates in LBC model
}
\examples{
iZZ <- matrix(0,100,10)
ixx <- matrix(0,100,1)
iyy <- matrix(0,100,1)
iw <-matrix(0,100,1)
# ixx: continuous predictor
# iyy: binary outcome
# iZZ: covariates to be incorporated in the model
# iw:  The weighted parameter
# myLBC <-
#   maxLik(
#     logLik = LogLikeFun,
#     grad = ScoreFun,
#     start = inits,
#    ixx = ixx,
#     iyy = iyy,
#    iw = iw,
#   iZZ = as.matrix(iZZ)
#  )
# as.matrix returns all values of iZZ as a matrix
}
\author{
li Xing
}

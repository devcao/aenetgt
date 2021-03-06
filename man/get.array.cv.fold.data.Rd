% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aenetgt.R
\name{get.array.cv.fold.data}
\alias{get.array.cv.fold.data}
\title{Splits array testing data into crossvalidation data sets.}
\usage{
get.array.cv.fold.data(X, Y, Z, K)
}
\arguments{
\item{X}{Design matrix with first column a column of 1s.}

\item{Y}{Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.}

\item{Z}{Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.}

\item{K}{The number of crossvalidation folds; \code{K} may be changed internally if the number of arrays is not divisible by \code{K}.}
}
\value{
List of length \code{K} of the crossvalidation training and testing data sets.
}
\description{
Splits array testing data into crossvalidation data sets.
}
\examples{
# generate individual covariate values and disease statuses
N <- 160
data <- model1(N)
X <- data$X
Y.true <- data$Yi
Se <- c(.92,.95) # set row/col and individual assay sensitivity
Sp <- c(.97,.98) # set row/col and individual assay specificity
cj <- 4 # set size of master pools
# subject individuals to array testing
assay.data <- array.assay.gen(Y.true,Se,Sp,cj)
Z <- assay.data$Z
Y <- assay.data$Y
array.cv.fold.data <- get.array.cv.fold.data(X,Y,Z,K=5)
}

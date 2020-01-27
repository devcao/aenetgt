\name{aenetgt-package}
\alias{aenetgt-package}
\alias{aenetgt}
\docType{package}
\title{
\packageTitle{aenetgt}
}
\description{
\packageDescription{aenetgt}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{aenetgt}
\packageIndices{aenetgt}

This package contains functions required to compute the adaptive elastic net logistic regression estimator on group testing data as described in Gregory, Wang, and McMahan (2018).  The principal functions are \code{enetgt}, which computes the elastic net estimator with a weighted l1 norm, and \code{mlegt}, which computes the maximum likelihood estimator.  Computing the adaptive elastic net estimator requires obtaining an initial estimate, either from \code{enetgt} or \code{mlegt}, and using the initial estimator to specify the weights in a second fitting with the \code{enetgt} function. The \code{enetgt.grid} function fits the elastic net estimator with weighted l1 norm over a grid of tuning parameter values.  The \code{aic.bic.eric.enetgt.grid} function takes the output of the \code{enetgt.grid} function and chooses tuning parameter values according to the AIC, BIC, and ERIC criteria.  The functions \code{get.individual.cv.fold.data}, \code{get.masterpool.cv.fold.data}, \code{get.dorfman.cv.fold.data}, and \code{get.assay.cv.fold.data} split group testing data into crossvalidation training and testing folds.  The function \code{cv.enetgt.grid} takes the crossvalidation data folds returned by any of these four functions as well as output from \code{enetgt.grid} and chooses the tuning parameter values according to a likelihood-based crossvalidation procedure. The functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, and \code{array.assay.gen} are provided for generating group testing data under individual testing, master pool testing, Dorfman testing, and array testing, respectively.

}
\author{
\packageAuthor{aenetgt}

Maintainer: \packageMaintainer{aenetgt}
}
\references{
	
	 Gregory, K. B., Wang, D., & McMahan, C. S. (2019). Adaptive elastic net for group testing. \emph{Biometrics}.
	 
	Hui, F.K., Warton, D.I., & Foster, S.D. (2015). Tuning parameter selection for the adaptive lasso using ERIC. 
		\emph{Journal of the Americal Statistical Association}, \emph{110}(509), 262-269.

	McMahan, C.S., Tebbs, J.M., Hanson, T.E., & Bilder, C.R. (2017). Bayesian regression for group testing data. 
		\emph{Biometrics}.
}

\keyword{package}
\examples{

# generate covariate values and disease statuses for 300 individuals from model0:
data <- model0(300)
X <- data$X
Y.true <- data$Y
b.true <- data$b
b.true

# subject individuals to individual testing
Se <- c(.94) # individual testing sensitivities
Sp <- c(.95) # individual testing specificities
pt.data <- individual.assay.gen(Y.true,Se,Sp,cj=1)
Z <- pt.data$Z
Y <- pt.data$Y

# compute the mle on the group testing data:
mlegt.out <- mlegt(X,Y,Z,Se,Sp,delta=.01)
b.mle <- mlegt.out$b.mle
b.mle

# compute adaptive elastic net estimator over a grid of tuning parameter values
n.lambda <- 8
n.theta <- 2
enetgt.grid.out <- enetgt.grid(X,Y,Z,Se,Sp,n.lambda,n.theta,weights = 1/abs(b.mle[-1]),delta=.01)

# make choices of the tuning parameters according to the aic, bic, and eric criteria
aic.bic.eric.enetgt.grid.out <- aic.bic.eric.enetgt.grid(enetgt.grid.out,Z,X,Se,Sp,"dorfman")
b.aenet.aic <- aic.bic.eric.enetgt.grid.out$b.enet.aic
b.aenet.bic <- aic.bic.eric.enetgt.grid.out$b.enet.bic
b.aenet.eric <- aic.bic.eric.enetgt.grid.out$b.enet.eric
b.aenet.aic 
b.aenet.bic
b.aenet.eric

# make a choice of the tuning parameter using 3-fold crossvalidation:
cv.fold.data <- get.individual.cv.fold.data(X,Y,Z,K=3)
cv.enetgt.grid.out <- cv.enetgt.grid(cv.fold.data,"individual",Se,Sp,enetgt.grid.out$lambda.seq,
						enetgt.grid.out$theta.seq,weights=1/abs(b.mle[-1]),
						B.INIT=enetgt.grid.out$B.ENET,delta=.01)
b.aenet.cv <- enetgt.grid.out$B.ENET[cv.enetgt.grid.out$cv.ind[1],,cv.enetgt.grid.out$cv.ind[2]]
b.aenet.cv
}

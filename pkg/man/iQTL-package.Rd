\name{iQTL-package}
\alias{iQTL-package}
%\alias{iQTL}
\docType{package}
\title{
Genome-wide QTL Mapping via H-likelihood
}
\description{
This package provides an R implementation of the double hierarchical generalized linear model (DHGLM) particularly for analyzing genome-wide marker data.
}
\details{
\tabular{ll}{
Package: \tab iQTL\cr
Type: \tab Package\cr
Version: \tab 2010.1\cr
Date: \tab 2010-07-26\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Depends: \tab VGAM\cr
}
% ~~ An overview of how to use the package, including the most important ~~
% ~~ functions ~~
}
\author{
Xia Shen

Maintainer: Xia Shen <xia.shen@lcb.uu.se>
}
\references{
Xia Shen, Lars Ronnegard and Orjan Carlborg (2010): Hierarchical likelihood opens a new way of estimating genetic values using dense marker maps. (14th QTLMAS, Poznan, Poland) SUBMITTED.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{GWAS, QTL, DHGLM, h-likelihood}
\seealso{
\code{\link{h.GWAS}}, \code{\link{h.testQTL}}
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
\dontrun{

## data

data(snpid)
data(pt)
data(gt)
data(pedi)
data(info)
marker.id <- info[snpid,1]

## set parameters

rho <- .9

## create spatial correlation matrix
## note that there are 5 chromosomes

border.idx <- c(0, 1971, 4032, 6080, 8068, 10029)
nmarker.chr <- numeric(5)
for (i in 1:5) nmarker.chr[i] <- sum((marker.id > border.idx[i])*
	(marker.id <= border.idx[i + 1]))
ac.mat <- matrix(0, dim(gt)[2], dim(gt)[2])
blockborder <- c(0, cumsum(nmarker.chr))
for (i in 1:5) ac.mat[(blockborder[i] + 1):blockborder[i + 1],
	(blockborder[i] + 1):blockborder[i + 1]] <- 
	rho**(toeplitz(1:nmarker.chr[i]) - 1)

## fixed effect design matrix

X <- model.matrix(~factor(pedi[,4]))
dimnames(X) <- dimnames(pt) <- dimnames(gt) <- c(NULL, NULL)

## analyze QT

GWAS <- h.GWAS(y = pt[,2], X = X, Z = gt, rho = rho, ac.mat = ac.mat, 
	phi.start = 2, plotting = FALSE)

par(mfrow = c(2, 1))
plot(marker.id, GWAS$v, xlab = 'Marker index', ylab = 'Marker effect', 
	type = 'p', cex = .6, col = 4, bty = 'n')
abline(v = border.idx, lty = 2, col = 'green4')
plot(marker.id, GWAS$lambda, xlab = 'Marker index', 
	ylab = 'Marker dispersion', type = 'l', col = 2, bty = 'n')
abline(v = border.idx, lty = 2, col = 'green4')

}
}

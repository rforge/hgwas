\name{h.GWAS}
\alias{h.GWAS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genome-wide Association via H-likelihood
}
\description{
This is an R implementation of the double hierarchical generalized linear model (DHGLM) particularly for analyzing genome-wide marker data.
}
\usage{
h.GWAS(y, X, Z, family = 'normal', link = 'logit', X.disp = NULL, 
       ac.mat = NULL, rho = 0, phi.start = 1, lambda.start = 0.1, 
       alpha.start = 1, conv.crit = 1e-05, max.iter = 200, plotting = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a numeric vector of the phenotypes.}
  \item{X}{a design matrix of the fixed effects.}
  \item{Z}{an incidence matrix of the markers across the genome. Each element in the matrix gives the genotype of the corresponding individual at the corresponding marker. See 'Details' for more about the coding.}
  \item{family}{a string indicating the distribution family of the phenotypes. \code{'normal'} and \code{'binary'} is available. See 'Details' for more information.}
  \item{link}{a string giving the link function for binary phenotypes. \code{'logit'}, \code{'probit'}, and \code{'cloglog'} are implemented and only avalable when \code{family = 'binary'}.}
  \item{X.disp}{a design matrix of the fixed effects of the dispersion part of the model, i.e. for modeling the marker-specific variance.}
  \item{ac.mat}{an autocorrelation matrix for the markers, which is used for smoothing the marker specific variances. If \code{NULL}, an spatial correlation matrix is created with the parameter \code{rho}. See 'Details' for more information.}
  \item{rho}{a numeric value in [0, 1) defining the spatial correlation for two adjacent markers. Only activated when \code{ac.mat = NULL}. See Details for more information.}
  \item{phi.start}{a numeric value giving the IWLS starting value of the dispersion parameter of the phenotypes. When \code{family = 'binary'}, this parameter is always fixed to be 1.}
  \item{lambda.start}{a numeric value giving the IWLS starting value of the dispersion parameter of the marker effects.}
  \item{alpha.start}{a numeric value giving the IWLS starting value of the dispersion parameter of the marker-specific variances. See 'References'.}
  \item{conv.crit}{a numeric value giving the convergence tolerance of the IWLS algorithm.}
  \item{max.iter}{an integer restricting the maximum number of IWLS interations.}
  \item{plotting}{a logic value specifying whether a figure will be produced after each IWLS iteration. The figures record the convergence behaviors of marker effects and marker-specific variances and will be saved in a single PDF file.}
}
\details{
Typically, for an intercross, at a particular marker, there are three kinds of genotypes: AA, AB, and BB. Coding of \code{Z} can be simply set to be the number of 'A' copies, i.e. 2 for 'AA', 1 for 'AB', and 0 for 'BB'. For a backcross, 1 and 0 can be used for coding heterozygotes and homozygotes, respectively.

When \code{family = 'binary'}, a link function is applied for linearizing the binary phenotypes to a continuous scale. The algorithm is not easy to converge depending on the starting values. \code{plotting = TRUE} is recommended to observe results for a certain number of iterations.

It is recommended that \code{ac.mat} is created according to the linkage disequilibium for each pair of markers across the genome. If \code{ac.mat = NULL}, the algorithm will create a spatial correlation matrix using argument \code{rho}. By default, \code{rho = 0} indicating no autocorrelation for the marker-specific variances. See 'References'.
}
\value{
  \item{phi}{the dispersion parameter of the phenotypes. When \code{family = 'binary'}, this parameter is always fixed to be 1.}
  \item{alpha}{the dispersion parameter of the marker-specific variances. See 'References'.}
  \item{lambda}{the marker-specific variances.}
  \item{beta}{the estimated fixed effects of the phenotypes.}
  \item{gamma}{the estimated fixed effects of the marker-specific variances.}
  \item{v}{the marker effects.}
  \item{b}{the estimated random effects of marker-specific variances.}
  \item{niter}{the number of iterations taken in the IWLS algorithm.}
}
\references{
Xia Shen, Lars Ronnegard and Orjan Carlborg (2010): Hierarchical likelihood opens a new way of estimating genetic values using dense marker maps. (14th QTLMAS, Poznan, Poland) SUBMITTED.
}
\author{
Xia Shen
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{iQTL-package}}
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


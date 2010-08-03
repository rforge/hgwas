\name{h.testQTL}
\alias{h.testQTL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Score test for QTL existence
}
\description{
This function performs QTL existence test via a score statistic.
}
\usage{
h.testQTL(y, Z, sig.level = .05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a numeric vector of the phenotypes.}
  \item{Z}{an incidence matrix of the markers across the genome. Each element in the matrix gives the genotype of the corresponding individual at the corresponding marker. See 'Details' for more about the coding.}
  \item{sig.level}{a numeric value between 0 and 1 specifying the significance level of QTL existence test.}
}
\details{
Typically, for an intercross, at a particular marker, there are three kinds of genotypes: AA, AB, and BB. Coding of \code{Z} can be simply set to be the number of 'A' copies, i.e. 2 for 'AA', 1 for 'AB', and 0 for 'BB'. For a backcross, 1 and 0 can be used for coding heterozygotes and homozygotes, respectively.
}
\value{
Screen print-out giving the score test statistic value and the p-value.
}
%\references{
%}
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

data(pt)
data(gt)

## test QTL existence using the last 10 F3 individuals

h.testQTL(pt[2317:2326,2], gt[2317:2326,]) # p-value = 0 < 0.05

## change the phenotypes and re-test

set.seed(123)
h.testQTL(rnorm(10), gt[2317:2326,]) # p-value = 0.7383826 > 0.05

}
}


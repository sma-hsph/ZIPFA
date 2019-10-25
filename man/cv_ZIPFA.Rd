\name{cv_ZIPFA}
\alias{cv_ZIPFA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross validation for Zero Inflated Poisson factor analysis
}
\description{
To conduct a cross validation for Zero Inflated Poisson factor analysis to find the number of factors.
}
\usage{
cv_ZIPFA(X, k, fold, tau = 0.1, cut = 0.8, tolLnlikelihood = 5e-4,
          iter = 20, tol = 1e-4, maxiter = 100, initialtau = 'iteration',
          Madj = TRUE, display = TRUE, parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{X}{The matrix to be decomposed.}
   \item{k}{A vector containing the number of factors to try.}
   \item{fold}{The number of folds used in cross validation.}
  \item{tau}{Initial tau value to fit. Will be overwritten by the first value in \code{initial} argument.}
  \item{cut}{To delete columns that has more than 100('\code{Cut}')\% zeros. \code{Cut = 1}, if no filtering.}
  \item{tolLnlikelihood}{The max percentage of log likelihood differences in two iterations.}
  \item{iter}{Max iterations.}
  \item{initialtau}{A character specifying the way to choose the initial value of tau at the beginning of EM iteration. \code{stable}: estimate tau from fitted beta in last round; \code{initial}: always use the initially assigned tau in \code{tau} or \code{initial}. Use the default tau = 0.1 if 'initial' is empty. \code{iteration}: use fitted tau in last round.}
  \item{tol}{Percentage of l2 norm change of [tau beta].}
  \item{maxiter}{Max iteration number in the zero inflated poisson regression.}
  \item{Madj}{If TRUE then adjust for relative library size M.}
  \item{display}{If TRUE display the fitting procedure.}
  \item{parallel}{Use \code{doParallel} and \code{foreach} package to accelerate.}
}
\details{
The function conducts cross validation on the zero-inflated Poisson factor analysis to determine the rank.
}
\value{
The function returns a matrix. Each row the CV likelihood of one fold. Each column is the result of number of factors in \code{k}.
}

\author{
Tianchen Xu
}


\examples{

\donttest{
data(X)
cv_ZIPFA(X, fold = 10, k = c(3,4))
}

}

\keyword{cross validation}


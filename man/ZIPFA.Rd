\name{ZIPFA}
\alias{ZIPFA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Zero Inflated Poisson factor analysis
}
\description{
To conduct a Zero Inflated Poisson factor analysis.
}
\usage{
ZIPFA(X, k, tau = 0.1, cut = 0.8, tolLnlikelihood = 5e-4,
        iter = 20, tol = 1e-4, maxiter = 100, initialtau = 'iteration',
        Madj = TRUE, display = TRUE, missing = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{X}{The matrix to be decomposed.}
   \item{k}{The number of factors.}
  \item{tau}{Initial tau value to fit. Will be overwritten by the first value in \code{initial} argument.}
  \item{cut}{To delete columns that has more than 100('\code{Cut}')\% zeros. \code{Cut = 1}, if no filtering.}
  \item{tolLnlikelihood}{The max percentage of log likelihood differences in two iterations.}
  \item{iter}{Max iterations.}
  \item{initialtau}{A character specifying the way to choose the initial value of tau at the beginning of EM iteration. \code{stable}: estimate tau from fitted beta in last round; \code{initial}: always use the initially assigned tau in \code{tau} or \code{initial}. Use the default tau = 0.1 if 'initial' is empty. \code{iteration}: use fitted tau in last round.}
  \item{tol}{Percentage of l2 norm change of [tau beta].}
  \item{maxiter}{Max iteration number in the zero inflated poisson regression.}
  \item{Madj}{If TRUE then adjust for relative library size M.}
  \item{display}{If TRUE display the fitting procedure.}
    \item{missing}{Reserved for \code{cv_ZIPFA}.}
}
\details{
The function conducts a zero-inflated Poisson factor analysis where the underlying Poisson rate is negatively associated with true zero probability.
}
\value{
\item{tau}{Fitted tau value.}
\item{Ufit}{A list containing fitted U matrix in each iteration. The last one is the final fit. }
\item{Vfit}{A list containing fitted V matrix in each iteration. The last one is the final fit. }
\item{itr}{Number of iterations.}
\item{Likelihood}{The likelihood for the training data.}
\item{CVLikelihood}{The likelihood for the testing data (if applicable)}
}

\author{
Tianchen Xu
}

\examples{

\donttest{
data(X)
ZIPFA(X, k = 3)}

}

\keyword{factor analysis}
\keyword{zero inflated}

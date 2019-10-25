\name{EMzeropoisson_mat}
\alias{EMzeropoisson_mat}

\title{Zero Inflated Possion Regression}
\description{
The zero inflated possion regression model.
}
\usage{
EMzeropoisson_mat(data, tau = 0.1, initial = NULL, initialtau = 'iteration',
                  tol = 1e-4, maxiter = 100, Madj = FALSE, m = NULL,
                  display = TRUE, intercept = TRUE)
}

\arguments{
  \item{data}{A matrix with the first columns is y and the rest columns are x.}
  \item{tau}{Initial tau value to fit. Will be overwritten by the first value in \code{initial} argument.}
  \item{initial}{A list of initial values for the fitting. \code{c(tau beta)}.}
  \item{initialtau}{A character specifying the way to choose the initial value of tau at the beginning of EM iteration. \code{stable}: estimate tau from fitted beta in last round; \code{initial}: always use the initially assigned tau in \code{tau} or \code{initial}. Use the default tau = 0.1 if 'initial' is empty. \code{iteration}: use fitted tau in last round.}
  \item{tol}{Percentage of l2 norm change of [tau beta].}
  \item{maxiter}{Max iteration number.}
  \item{Madj}{If TRUE then adjust for relative library size M.}
  \item{m}{A vector containing relative library size M.}
  \item{display}{If TRUE display the fitting procedure.}
  \item{intercept}{If TRUE then the model contains an intercept.}
}
\details{
The function estimates the coefficients in a new type of zero-inflated Poisson regression where the underlying Poisson rate is negatively associated with true zero probability.
}
\value{
The function turns a matrix. Each row is fitted value in each iteration. The last row the final result. The first column is fitted tau. If \code{intercept} is ture, then the second column is the intercept, and the rest columns are other coefficients. If \code{intercept} is false, the rest columns are other coefficients.
}

\author{
Tianchen Xu
}

\examples{
n = 5000;
x1 = rnorm(n);
x2 = rnorm(n);
lam = exp(x1 - 2*x2 + 1.5);
y = rpois(n, lam)
tau = .75
p = 1/(1+lam^tau);
Z = rbinom(n, 1, p);
y[as.logical(Z)] = 0;

res = EMzeropoisson_mat(matrix(c(y,x1,x2),ncol=3), Madj = FALSE, intercept = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{factor analysis}% use one of  RShowDoc("KEYWORDS")
\keyword{zero inflated}

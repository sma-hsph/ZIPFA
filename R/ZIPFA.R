#' @import Matrix
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import optimx
#' @import trustOptim
#' @importFrom stats glm.fit poisson median dpois
#' @importFrom utils capture.output
#' @export cv_ZIPFA
#' @export ZIPFA
#' @export EMzeropoisson_mat


cv_ZIPFA <- function(X, k, fold, tau = 0.1, cut = 0.8, tolLnlikelihood = 5e-4, iter = 20, tol = 1e-4, maxiter = 100, initialtau = 'iteration', Madj = TRUE, display = TRUE, parallel = FALSE){
  # cv_ZIPFA    To conduct cross validation on ZIPFA model.
  #   [cvsample,Allres] = cv_ZIPFA (X, k)
  #
  #   [cvsample,Allres] = cv_ZIPFA (X, k, fold, tau)
  #
  #   [cvsample,Allres] = cv_ZIPFA (X, k, 'cut', 1)
  #
  #   [cvsample,Allres] = cv_ZIPFA (X, k, 'display', false, 'savemat', true, ...)
  #
  # -X: The matrix to be decomposed.
  # -k: The number of factors. It can be a vector.
  # -fold (10): The number of folds used in cross validation.
  # -tau (0.1): The initial guess for tau.
  # -'cut' (0.8): Whether to delete columns that has more than 100('cut')% zeros. 'Cut' = 1, if no filtering.
  # -'tolLnlikelihood' (5e-4): The max percentage of log likelihood differences in two iterations.
  # -'iter' (20): Max iterations.
  # -'tol' (1e-4): Percentage of l2 norm change of [tau beta] in ZIP regression.
  # -'maxiter' (100): Max iterations in ZIP regression.
  # -'initialtau' (iteration'): Choose the initial value of tau at the beginning of EM iteration in ZIP regression.
  #     'stable': estimate tau from fitted beta in last round;
  #     'initial': always use the initially assigned tau in 'tau' or 'initial';
  #           Use the default tau = 0.1 if 'initial' is empty.
  #     'iteration': use fitted tau in last round.
  # -'Madj' (true): Whether adjust for relative library size M.
  # -'display' (true): Display the fitting procedure. Info in ZIPFA will not be shown in 'Parallel' mode even 'Display' is true.
  # -'savemat' (false): Whether to save ZIPFA results in all factor numbers and each fold.
  # -'parallel' (true): Use parallel toolbox to accelerate.


  if (!is.matrix(X)) stop('data should be a matrix.')
  if (any(is.na(X))) stop('There is Nas in X !')
  m <- nrow(X)
  n0 <- ncol(X)
  if (!(is.atomic(tau) && length(tau) == 1L)) stop('tau should a scalar.')
  if (!any(initialtau == c('stable','initial','iteration'))) stop('initialtau must be one of "statble", "initial", "iteration".')

  X <- X[, colSums(X==0)/m < cut]
  n <- ncol(X)
  if (n!=n0){
    warning(sprintf('%d columns have been dropped due to the "cut" option.\n',n0-n))
  }

  # generate assignment
  cvsample <- NULL
  for (i in 1:(fold-1)){
    cvsample <- c(cvsample, i * rep(1, floor(m*n/fold)))
  }
  cvsample <- c(cvsample, fold * rep(1, m*n - (fold-1) * floor(m*n/fold)))
  cvsample <- sample(cvsample, length(cvsample), replace = FALSE)
  cvsample <- matrix(cvsample, nrow = m, ncol = n)

  if (exists('cvsample' , envir = .GlobalEnv, inherits = FALSE)){
    message('Use existing cvsample in the global environment.\n')
    cvsample <- get('cvsample', envir = .GlobalEnv, inherits = FALSE)
  }
  # else{
  #   assign('cvsample', cvsample, envir = .GlobalEnv, inherits = FALSE)
  # }
  if (parallel){
    cl <- parallel::makeCluster(detectCores(logical = FALSE))
    doParallel::registerDoParallel(cl)
  }

  Allres <- NULL
  for (nf in k){
    Mlike <- matrix(0, nrow = fold, ncol = 1)
    if (parallel){
      Mres <- foreach(rept=1:fold, .combine = 'c', .packages = c('Matrix', 'trustOptim', 'optimx')) %dopar% {
        missing <- cvsample
        missing[missing==rept] <- 0
        missing[missing!=0] <- 1
        resfit <- ZIPFA(X = X, k = nf, tau = 0.1, cut = cut, tolLnlikelihood = tolLnlikelihood, iter = iter, tol = tol, maxiter = maxiter, initialtau = initialtau, Madj = Madj, display = display, missing = (missing==1))
        if (display){
          message(sprintf('\n factor = %d; fold = %d; esttau = %.3g; Likelihood = %.5g; MLikelihood = %.4g \n',nf,rept,resfit$tau,resfit$Likelihood[length(resfit$Likelihood)],resfit$CVLikelihood[length(resfit$CVLikelihood)]))
        }
        resfit$CVLikelihood[length(resfit$CVLikelihood)]
      }
      Mlike[1:fold,1] <- Mres
    }else{
      for (rept in 1:fold){
        missing <- cvsample
        missing[missing==rept] <- 0
        missing[missing!=0] <- 1
        resfit <- ZIPFA(X = X, k = nf, tau = 0.1, cut = cut, tolLnlikelihood = tolLnlikelihood, iter = iter, tol = tol, maxiter = maxiter, initialtau = initialtau, Madj = Madj, display = display, missing = (missing==1))
        if (display){
          message(sprintf('\n factor = %d; fold = %d; esttau = %.3g; Likelihood = %.5g; MLikelihood = %.4g \n',nf,rept,resfit$tau,resfit$Likelihood[length(resfit$Likelihood)],resfit$CVLikelihood[length(resfit$CVLikelihood)]))
        }
        Mlike[rept, 1] <- resfit$CVLikelihood[length(resfit$CVLikelihood)]
      }
    }
    colnames(Mlike) <- paste0('K_', nf)
    Allres <- cbind(Allres, Mlike)
  }
  parallel::stopCluster(cl = cl)
  return(Allres)
}


ZIPFA <- function(X, k, tau = 0.1, cut = 0.8, tolLnlikelihood = 5e-4, iter = 20, tol = 1e-4, maxiter = 100, initialtau = 'iteration', Madj = TRUE, display = TRUE, missing = NULL){
  # ZIPFA   To conduct the zero inflated zero factor analysis.
  #   [Ufit,Vfit] = ZIPFA(X, k)
  #
  #   [Ufit,Vfit] = ZIPFA(X, k, tau, 'cut', 1, 'display', false, ...)
  #
  #   [Ufit,Vfit,itr,finaltau,Likelihood] = ZIPFA(X, k)
  #
  #   [Ufit,Vfit,itr,finaltau,Likelihood,CVLikelihood] = ZIPFA(X, k, 'missing', missingmat)
  #
  # -X: The matrix to be decomposed.
  # -k: The number of factors.
  # -tau (0.1): The initial guess for tau.
  # -'cut' (0.8): Whether to delete columns that has more than 100('cut')% zeros. 'Cut' = 1, if no filtering.
  # -'tolLnlikelihood' (5e-4): The max percentage of log likelihood differences in two iterations.
  # -'iter' (20): Max iterations.
  # -'tol' (1e-4): Percentage of l2 norm change of [tau beta] in ZIP regression.
  # -'maxiter' (100): Max iterations in ZIP regression.
  # -'initialtau' (iteration'): Choose the initial value of tau at the beginning of EM iteration in ZIP regression.
  #     'stable': estimate tau from fitted beta in last round;
  #     'initial': always use the initially assigned tau in 'tau' or 'initial';
  #           Use the default tau = 0.1 if 'initial' is empty.
  #     'iteration': use fitted tau in last round.
  # -'Madj' (true): Whether adjust for relative library size M.
  # -'display' (true): Display the fitting procedure.
  # -'missing' ([]): TRUE/FALSE matrix. If 'missing' is not empty, then CVLikelihood is likelihood of X with missing = TRUE.

  if (!is.matrix(X)) stop('data should be a matrix.')
  m <- nrow(X)
  n0 <- ncol(X)
  if (!(is.atomic(tau) && length(tau) == 1L)) stop('tau should a scalar.')
  if (!(is.atomic(k) && length(k) == 1L)) stop('k should a scalar.')
  if (!any(initialtau == c('stable','initial','iteration'))) stop('initialtau must be one of "statble", "initial", "iteration".')

  if(!is.null(missing)){
    missing <- missing[, colSums(X==0)/m < cut]
  }
  X <- X[, colSums(X==0)/m < cut]
  n <- ncol(X)
  if (n!=n0){
    warning(sprintf('%d columns have been dropped due to the "cut" option.\n',n0-n))
  }

  Xt <- t(X)
  if (any(is.na(X))){
    if (!is.null(missing)){
      stop('Cannot assign "missing" when X already has missing values')
    }
    indexram <- !is.na(as.vector(X))
    indexram_row <- !is.na(as.vector(Xt))
  }else{
    if (!is.null(missing)){
      indexram <- as.vector(missing)
      indexram_row <- as.vector(t(missing))
    }else{
      indexram <- indexram_row <- 1:(m*n)
    }
  }

  XNA <- X
  XNA[XNA==0] <- NA
  lambdahat <- colMeans(log(XNA), na.rm = TRUE)
  for (i in 1:n){
    XNA[is.na(XNA[,i]), i] <- round(exp(lambdahat[i]))
  }

  presvd <- svd(log(XNA))
  Uold <- presvd$u[,1:k,drop=FALSE] %*% diag(presvd$d[1:k])
  Vold <- presvd$v[,1:k,drop=FALSE]

  Ufit <- vector('list', iter)
  Vfit <- vector('list', iter)
  Likelihood <- NULL
  CVLikelihood <- NULL

  if (Madj){
    ratio <- rep(rowSums(X, na.rm = TRUE)/median(rowSums(X, na.rm = TRUE)), n)
    ratio_row <- rep(rowSums(X, na.rm = TRUE)/median(rowSums(X, na.rm = TRUE)), each = n)
  }

  for (itr in 1:iter){
    if (display){
      message(sprintf('\n ****************** \n Round %.0f \n ******************\n', itr))
    }

    # updata U
    if (Madj){
      mi <- ratio_row[indexram_row]
    }else{
      mi <- 1
    }
    Voldc <- bdiag(rep(list(Vold), m))
    dat <- cbind(matrix(as.vector(Xt), ncol=1), Voldc)
    dat <- dat[indexram_row, , drop=FALSE]
    if (itr==1){
      uuuu <- EMzeropoisson_mat(data = dat, tau = 0.1, initial = cbind(tau, matrix(as.vector(t(Uold)), nrow = 1)), initialtau = initialtau, tol = tol, maxiter = maxiter, Madj = Madj, m = matrix(mi, ncol = 1), display = display, intercept = FALSE)
    }else{
      uuuu <- EMzeropoisson_mat(data = dat, tau = 0.1, initial = cbind(uuuu[nrow(uuuu),1], matrix(as.vector(t(Ufit[[itr-1]])), nrow = 1)), initialtau = initialtau, tol = tol, maxiter = maxiter, Madj = Madj, m = matrix(mi, ncol = 1), display = display, intercept = FALSE)
    }

    Unew <- matrix(uuuu[nrow(uuuu),-1], byrow = TRUE, ncol = k)

    # update V
    if (Madj){
      mi <- ratio[indexram]
    }else{
      mi <- 1
    }
    Uoldc <- bdiag(rep(list(Uold), n))
    dat <- cbind(matrix(as.vector(X), ncol=1), Uoldc)
    dat <- dat[indexram, , drop=FALSE]
    if (itr==1){
      uuuu <- EMzeropoisson_mat(data = dat, tau = 0.1, initial = cbind(tau, matrix(as.vector(t(Vold)), nrow = 1)), initialtau = initialtau, tol = tol, maxiter = maxiter, Madj = Madj, m = matrix(mi, ncol = 1), display = display, intercept = FALSE)
    }else{
      uuuu <- EMzeropoisson_mat(data = dat, tau = 0.1, initial = cbind(uuuu[nrow(uuuu),1], matrix(as.vector(t(Vfit[[itr-1]])), nrow = 1)), initialtau = initialtau, tol = tol, maxiter = maxiter, Madj = Madj, m = matrix(mi, ncol = 1), display = display, intercept = FALSE)
    }

    Vnew <- matrix(uuuu[nrow(uuuu),-1], byrow = TRUE, ncol = k)

    if (Madj){
      Likelihood <- c(Likelihood, suppressMessages(lnL_mat(data = cbind(matrix(as.vector(X), ncol=1), Uoldc), coefficient = uuuu, subset = indexram, mi = ratio)))
      if (!is.null(missing)){
        CVLikelihood <- c(CVLikelihood, suppressMessages(lnL_mat(data = cbind(matrix(as.vector(X), ncol=1), Uoldc), coefficient = uuuu, subset = !indexram, mi = ratio)))
      }
    }else{
      Likelihood <- c(Likelihood, suppressMessages(lnL_mat(data = cbind(matrix(as.vector(X), ncol=1), Uoldc), coefficient = uuuu, subset = indexram, mi = NULL)))
      if (!is.null(missing)){
        CVLikelihood <- c(CVLikelihood, suppressMessages(lnL_mat(data = cbind(matrix(as.vector(X), ncol=1), Uoldc), coefficient = uuuu, subset = !indexram, mi = NULL)))
      }
    }

    if (display){
      message(sprintf('\n Likelihood is %f \n',Likelihood))
      if (!is.null(CVLikelihood)){
        message(sprintf('CV Likelihood is %f \n',CVLikelihood))
      }
    }

    # next step
    presvd <- svd(Unew %*% t(Vnew))
    Uold <- presvd$u[,1:k,drop=FALSE] %*% diag(presvd$d[1:k])
    Vold <- presvd$v[,1:k,drop=FALSE]

    # save answer
    Ufit[[itr]] <- Uold
    Vfit[[itr]] <- Vold

    if (itr!=1){
      if (display){
        message(sprintf('Max Ln likelihood diff = %.4g %% \n',100*(Likelihood[itr]-Likelihood[itr-1])/abs(Likelihood[itr-1])))
      }
      if (((Likelihood[itr]-Likelihood[itr-1])/abs(Likelihood[itr-1]))<tolLnlikelihood) break
    }
  }

  finaltau <- uuuu[nrow(uuuu), 1]
  Ufit <- Ufit[!sapply(Ufit, is.null)]
  Vfit <- Vfit[!sapply(Vfit, is.null)]
  return(list(tau = finaltau, Ufit = Ufit, Vfit = Vfit, itr = itr, Likelihood = Likelihood, CVLikelihood = CVLikelihood))
}


EMzeropoisson_mat <- function(data, tau = 0.1, initial = NULL, initialtau = 'iteration', tol = 1e-4, maxiter = 100, Madj = FALSE, m = NULL, display = TRUE, intercept = TRUE){
  # EMzeropoisson_mat  To fit the zero inflated Poisson regression.
  #   fittedbeta = EMzeropoisson_mat([y x])
  #
  #   fittedbeta = EMzeropoisson_mat([y x], tau, 'display', false, ...)
  #
  # -data: First y then x.
  # -tau (0.1): Initial tau to fit. Will be overwritten by 'initial'.
  # -'initial' ([]): Initial [tau beta].
  # -'initialtau' (iteration'): Choose the initial value of tau at the beginning of EM iteration.
  #     'stable': estimate tau from fitted beta in last round;
  #     'initial': always use the initially assigned tau in 'tau' or 'initial';
  #         Use the default tau = 0.1 if 'initial' is empty.
  #     'iteration': use fitted tau in last round.
  # -'tol' (1e-4): Percentage of l2 norm change of [tau beta].
  # -'maxiter' (100): Max iteration.
  # -'Madj' (false): Whether adjust for relative library size M.
  # -'m' ([]): Relative library size M.
  # -'display' (true): Display the fitting procedure.
  # -'intercept' (true): Whether the model contains an intercept.
  #
  # Result contains the fitted results in each row. The last row shows
  # the final result. First column is tau, second column is intercept (if the
  # model has intercept), other columns are fitted coefficients.

  n <- nrow(data)
  # if (!is.matrix(data)) stop('data should be a matrix.')
  if (!(is.atomic(tau) && length(tau) == 1L)) stop('tau should a scalar.')
  if (!any(initialtau == c('stable','initial','iteration'))) stop('initialtau must be one of "statble", "initial", "iteration".')
  if (Madj == TRUE && !(ncol(m)==1 & nrow(m)==n)) stop('m should match y.')

  x <- data[, -1]
  y <- data[, 1, drop=FALSE]



  if (intercept){
    Dx <- cbind(matrix(1, nrow = n, ncol = 1), x)
  }else{
    Dx <- x
  }

  if (is.null(initial)){
    if (display){
      message('Initialzing... \n')
    }
    startvar <- glm.fit(Dx, y, family = poisson(link = "log"))$coefficients
    startvar <- matrix(c(tau, startvar), nrow = 1)
  }else{
    startvar <- initial
  }

  i <- 1
  if (display) message('Start maximizing... \n')
  while ((i < maxiter) & (i==1 || (sqrt(sum((startvar[i,]-startvar[i-1,])^2))/sqrt(sum(startvar[i,]^2))>tol))){
    # E step
    Z <- matrix(0, nrow = 1, ncol = n)
    Izero <- as.vector(y==0)
    if (Madj){
      Z[1, Izero] <- as.vector(1/(1 + exp(startvar[i, 1] * (Dx[Izero, , drop=FALSE] %*% t(startvar[i, -1, drop=FALSE]))-exp(m[Izero, 1, drop=FALSE] * (Dx[Izero,, drop=FALSE] %*% t(startvar[i, -1, drop=FALSE]))) )))
    }else{
      Z[1, Izero] <- as.vector(1/(1 + exp(startvar[i, 1] * (Dx[Izero, , drop=FALSE] %*% t(startvar[i, -1, drop=FALSE]))-exp((Dx[Izero,, drop=FALSE] %*% t(startvar[i, -1, drop=FALSE]))) )))
    }

    # M step
    if (i!=1 & initialtau=='stable'){
      tau <- -log( n/ sum(Z==0) - 1) / mean((Dx %*% t(startvar[i, -1, drop=FALSE])))
    }else if (initialtau=='initial'){
      tau <- startvar[1, 1]
    }else if (initialtau=='iteration'){
      tau <- startvar[i, 1]
    }

    if (Madj){
      # optimx(c(tau, startvar[i,-1]), fn = value, gr = grad, hess = hess, method = c('BFGS'), control=list(save.failures=TRUE, trace=0, kkt = FALSE, starttests = FALSE, dowarn = FALSE, usenumDeriv = FALSE))
      nnnn <- capture.output(resfit <- trust.optim(c(tau, startvar[i,-1]), fn = value, hs = hess, gr = grad, method =  "Sparse", control = list(report.level=0), rmc = m, Dx=Dx, y=y, Z=Z))
      if (resfit$status %in% c('Success','Radius of trust region is less than stop.trust.radius')){
        startvar <- rbind(startvar, matrix(resfit$solution, nrow=1))
      }else{
        warning('trust.optim fails, try BFGS')
        resfit <- optimx(c(tau, startvar[i,-1]), fn = value, gr = grad, hess = hess, method = 'BFGS', control=list(save.failures=TRUE, trace=0, kkt = FALSE, starttests = FALSE, dowarn = FALSE, usenumDeriv = FALSE), rmc=m, Dx=Dx, y=y, Z=Z)
        startvar <- rbind(startvar, as.matrix(resfit[1, 1:ncol(startvar)]))
      }
    }else{
      nnnn <- capture.output(resfit <- trust.optim(c(tau, startvar[i,-1]), fn = value, hs = hess, gr = grad, method =  "Sparse", control = list(report.level=0), rmc = 1, Dx=Dx, y=y, Z=Z))
      if (resfit$status %in% c('Success','Radius of trust region is less than stop.trust.radius')){
        startvar <- rbind(startvar, matrix(resfit$solution, nrow=1))
      }else{
        warning('trust.optim fails, try BFGS')
        resfit <- optimx(c(tau, startvar[i,-1]), fn = value, gr = grad, hess = hess, method = 'BFGS', control=list(save.failures=TRUE, trace=0, kkt = FALSE, starttests = FALSE, dowarn = FALSE, usenumDeriv = FALSE), rmc=1, Dx=Dx, y=y, Z=Z)
        startvar <- rbind(startvar, as.matrix(resfit[1, 1:ncol(startvar)]))
      }
    }

    i <- i+1
    if (initialtau =='iteration' & max(abs(diff(startvar[,1])), na.rm = TRUE)>50 & i>3){
      warning('May be divergent tau. Try "stable=" model.')
    }
    if (display){
      message(sprintf('This is %.0f th iteration, Frobenius norm diff = %g. \n', i, sqrt(sum((startvar[i,]-startvar[i-1,])^2))/sqrt(sum(startvar[i,]^2))))
    }
  }
  return(startvar)
}




value <- function(beta, rmc = NULL, Dx, y, Z){
  tau <- beta[1]
  beta <- matrix(beta[-1], ncol = 1)
  xb <- Dx %*% beta
  mexb <- rmc * exp(xb)
  value <- tau * Z %*% xb + sum(log(1+exp(-tau * xb))) + (Z - 1) %*% (y * xb - mexb)
  value <- as.numeric(value)
  if (is.infinite(value) | is.na(value)){
    small <- (-tau * xb)<30
    value <- tau * Z %*% xb + sum(log(1+exp(-tau * xb)) * small + (-tau * xb) * (!small), na.rm = TRUE) + (-1 + Z) %*% (y * xb - mexb)
    value <- as.numeric(value)
  }
  return(value)
}

grad <- function(beta, rmc = NULL, Dx, y, Z){
  tau <- beta[1]
  beta <- matrix(beta[-1], ncol = 1)
  xb <- Dx %*% beta
  mexb <- rmc * exp(xb)
  etauxb <- exp(tau * xb)
  wz <- t(1/(etauxb + 1)) - Z
  dldb <- (wz * tau + (1 - Z) * t(y-mexb)) %*% Dx
  dldtau <- wz %*% xb;
  # return(-t(cbind(dldtau, dldb)))
  # return(-as.numeric(cbind(dldtau, dldb)))
  return(-c(as.numeric(dldtau), as.numeric(dldb)))
}

hess <- function(beta, rmc = NULL, Dx, y, Z){
  tau <- beta[1]
  beta <- matrix(beta[-1], ncol = 1)
  xb <- Dx %*% beta
  etauxb <- exp(tau * xb)
  dldtautau <- -sum(xb^2 * etauxb/(etauxb + 1)^2);
  dldtaudb <- (t((etauxb + 1 - tau * xb * etauxb)/(etauxb + 1)^2) - Z) %*% Dx
  tp <- tau * tau * t(-etauxb/(etauxb + 1)^2) + (1 - Z) * t(-rmc * exp(xb))
  dldbdb <- t(Dx * as.vector(tp)) %*% Dx
  hess <- -rbind(cbind(dldtautau, dldtaudb), cbind(t(dldtaudb), dldbdb))
  if (any(is.infinite(hess)) | any(is.na(hess))){
    small <- etauxb < exp(30)
    dldtautau <- -sum( xb^2 * etauxb/(etauxb + 1)^2 * small + xb^2/etauxb * (!small), na.rm = TRUE)
    v1 <- (etauxb + 1 - tau * xb * etauxb)/(etauxb + 1)^2
    v2 <- (1 - tau * xb)/etauxb
    V <- Matrix(0, nrow = nrow(y), ncol = 1)
    V[small] <- v1[small]
    V[!small] <- v2[!small]
    dldtaudb <- (t(V) - Z) %*% Dx
    R <- Matrix(0, nrow = nrow(y), ncol = 1)
    r1 <- etauxb/(etauxb+1)^2
    r2 <- 1/etauxb
    R[small] <- r1[small]
    R[!small] <- r2[!small]
    tp <- -tau * tau * t(R) + (1 - Z) * t(-rmc * exp(xb))
    dldbdb =  t(Dx * as.vector(tp)) %*% Dx
    hess <- -rbind(cbind(dldtautau, dldtaudb), cbind(t(dldtaudb), dldbdb))
  }
  # write.csv(hess,'hhh.csv')
  return(Matrix(hess, sparse = TRUE))
}

lnL_mat <- function(data, coefficient, subset, mi){
  coefficient <- coefficient[nrow(coefficient), , drop=FALSE]
  tau <- coefficient[1, 1]
  beta <- t(coefficient[1, -1, drop=FALSE])
  y <- data[, 1, drop=FALSE]
  n <- length(y)

  Dx <- data[, -1, drop=FALSE]
  lnlam <- Dx %*% beta
  p <- 1/(1 + exp(tau * lnlam))
  if (is.null(mi)){
    mi <- Matrix(1, nrow = n, ncol = 1)
  }else{
    mi <- Matrix(mi)
  }

  like <- y
  indexz <- y==0
  indexnotz <- y>=20
  indexnotzl <- ((y!=0) & (y<20))

  like[indexz] <- log(p[indexz]+(1-p[indexz])*(dpois(0,mi[indexz]*exp(lnlam[indexz]))))

  like[indexnotz] <- log(1-p[indexnotz]) - mi[indexnotz]*exp(lnlam[indexnotz]) + y[indexnotz]*(log(mi[indexnotz])+lnlam[indexnotz])-(0.5*log(2*pi*y[indexnotz]) + y[indexnotz]*(log(y[indexnotz])-1))

  like[indexnotzl] <- log(1-p[indexnotzl]) - mi[indexnotzl]*exp(lnlam[indexnotzl]) + y[indexnotzl]*(log(mi[indexnotzl])+lnlam[indexnotzl])-log(factorial(y[indexnotzl]))

  if (!is.null(subset)){
    like <- like[subset]
  }
  return(sum(like))
}


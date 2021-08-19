#' Estimate the CSTE curve for time to event outcome with right censoring.  
#'
#' Estimate the CSTE curve for time to event outcome with right censoring. 
#' The working model 
#' is \deqn{\lambda(t| X, Z) = \lambda_0(t) \exp(g_1(X\beta_1)Z + g_2(X\beta_2)),}
#' which implies that \eqn{CSTE(x) = g_1(x\beta_1)}.   
#' 
#'@param x samples of covariates which is a \eqn{n*p} matrix.
#'@param y samples of time to event which is a \eqn{n*1} vector.
#'@param z samples of treatment indicator which is a \eqn{n*1} vector.
#'@param status samples of censoring indicator which is a \eqn{n*1} vector, default value is NULL, indicating no censoring.  
#'@param beta_ini initial values for \eqn{(\beta_1', \beta_2')'}, default value is NULL.
#'@param lam value of the lasso penalty parameter \eqn{\lambda} for \eqn{\beta_1} and
#'\eqn{\beta_2}, default value is 0.
#'@param nknots number of knots for the B-spline for estimating \eqn{g_1} and \eqn{g_2}.
#'@param max.iter maximum iteration for the algorithm.
#'@param eps numeric scalar \eqn{\geq} 0, the tolerance for the estimation of \eqn{\beta_1} and \eqn{\beta_2}. 
#'
#'@return A S3 class of cste, which includes:
#' \itemize{
#'    \item \code{beta1}: estimate of \eqn{\beta_1}.
#'    \item \code{beta2}: estimate of \eqn{\beta_2}.
#'    \item \code{B1}: the B-spline basis for estimating \eqn{g_1}.
#'    \item \code{B2}: the B-spline basis for estimating \eqn{g_2}. 
#'    \item \code{delta1}: the coefficient of B-spline for estimating \eqn{g_1}.
#'    \item \code{delta2}: the coefficient for B-spline for estimating \eqn{g_2}.
#'    \item \code{iter}: number of iteration.
#'    \item \code{g1}: the estimate for \eqn{g_1(X\beta_1)}.
#'    \item \code{g2}: the estimate for \eqn{g_2(X\beta_2)}. 
#' }
#'@examples
#' ## Quick example for the cste
#' 
#' # --------  Example 1: Simulated Data ---------  # 
#' ## generate data 
#' set.seed(100)
#' X1 <- runif(200,0,20)
#' X2 <- runif(200,0,20)
#' X3 <- runif(200,0,20)
#' X <- cbind(X1, X2, X3)
#' Z <- c(rep(1,100),rep(0,100))
#' 
#' beta <- c(0.3,0.4,0.5)
#' lambda <- 0.1 * exp(log(X%*%beta/10+0.6)*Z+(X%*%beta)/20)
#' Time <- -log(runif(200,0,1))/lambda
#' C <- runif(200,12,15)
#' S <- as.numeric(Time <= C)
#' Time <- Time*S + C*(1-S)
#' 
#' ## estimate the CSTE curve
#' fit <- cste_surv(X, Time, Z, S, nknots=3, max.iter=200)
#' 
#' ## pointwise confidence interval (PCI)
#' res <- cste_surv_PCI(fit)
#' 
#' ## plot
#' plot(res$or_x, res$fit_x, col = 'red', type = "l", 
#'      xlim=c(8,18), ylim=c(-1,3), lwd = 2, 
#'      ylab = "CSTE", xlab = "X * beta1",
#'      main ="Pointwise Confidence interval")
#' lines(res$or_x, res$lower_bound, lwd = 3, col = 'purple', lty = 2)
#' lines(res$or_x, res$upper_bound, lwd = 3, col = 'purple', lty = 2)
#' abline(h = 0, lty = 2, cex = 0.2)
#' legend("topleft", legend = c("Estimates", "PCI"), 
#'         lwd = c(2,3), lty = c(1, 2), col = c('red','purple'))
#' # xb <- sort(as.numeric(X%*%beta))
#' # lines(xb, log(xb/10+0.6), type='l', col='blue')
#' 
#' 
#' # --------  Example 2: Real Data ---------  # 
#' dat <- survival::lung
#' 
#' ## estimate the CSTE curve
#' fit <- cste_surv(dat[,6:8], dat$time, dat$sex, 
#'                  dat$status, nknots = 3, max.iter = 100)
#'                  
#' ## 95% pointwise confidence interval (PCI)
#' res <- cste_surv_PCI(fit)
#' 
#' ## plot
#' plot(res$or_x, res$fit_x, col = 'red', type = "l", lwd = 2, ylim = c(-5,2),
#'      ylab = "CSTE", xlab = "X * beta1", main ="Pointwise Confidence interval")
#' lines(res$or_x, res$lower_bound, lwd = 3, col = 'purple', lty = 2)
#' lines(res$or_x, res$upper_bound, lwd = 3, col = 'purple', lty = 2)
#' abline(h = 0, lty = 2, cex = 0.2)
#' legend("topleft", legend = c("Estimates", "PCI"), 
#'         lwd = c(2,3), lty = c(1, 2), col = c('red','purple'))
#' 
#' @references
#' Ma Y. and Zhou X. (2017). 
#' Treatment selection in a randomized clinical trial via covariate-specific 
#' treatment effect curves, \emph{Statistical Methods in Medical Research}, 26(1), 124-141.
#' 
#' @seealso \code{\link{cste_surv_PCI}, \link{predict_cste_surv}, \link{select_cste_surv}}


# cste estimation for survival data 
cste_surv <- function(x, y, z, status = NULL, beta_ini = NULL, lam = 0, nknots = 2, max.iter = 200, eps = 1e-3) {
  cmplt <- complete.cases(x)
  x <- x[cmplt,]
  y <- y[cmplt]
  z <- z[cmplt]
  status <- status[cmplt]
  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  atrisk <- matrix(0,n,n)
  for (i in 1:n){
    atrisk[i,] <- as.numeric(y >= y[i])
  }
  if (is.null(status)) stutus<-rep(1,n)
  if(p==1) {
    B1 <- B2 <- bs(x, df = nknots+4, intercept = TRUE) 
    knots <- seq(0, 1, length = nknots + 2)
    B <- cbind(z * B1, B2)
    fit <- coxph(Surv(y, status) ~ 0+B)
    delta1 <- coef(fit)[1:(nknots+4)]
    delta2 <- coef(fit)[(nknots+5):(2*nknots+8)]
    delta1[is.na(delta1)] = 0
    delta2[is.na(delta2)] = 0
    geta <- z * B1 %*% delta1 + B2 %*% delta2
    g1 <- B1 %*% delta1
    g2 <- B2 %*% delta2
    loss <- -fit$loglik[2]
    bic <- -2 * loss + (nknots+4) * log(n) 
    aic <- -2 * loss + 2 * (nknots+4) 
    out <- list(beta1 = 1, beta2 = 1, B1 = B1, B2 = B2, delta2 = delta2, delta1 = delta1, g = geta, knots = knots, x = x, y = y, z = z, status = status, nknots = nknots, p=p, g1=g1, g2=g2, bic=bic, aic=aic)
    class(out) <- "cste"
    return(out)
  } else {
    flag <- FALSE
    # if(is.null(truth)) truth <- rep(p,2)
    truth <- rep(p,2)
    if(is.null(beta_ini)) beta_ini <- c(normalize(rep(1, truth[1])), normalize(rep(1, truth[2])))
    else beta_ini <- c(beta_ini[1:truth[1]], beta_ini[(truth[1]+1):sum(truth)])
    beta_curr <- beta_ini
    conv <- FALSE
    iter <- 0
    # knots location
    knots <- seq(0, 1, length = nknots + 2)
    len.delta <- length(knots) + 2
    while(conv == FALSE & iter < max.iter) {
      iter <- iter + 1
      # step 1 fix beta to estimate g
      beta1 <- beta_curr[1:truth[1]]
      beta2 <- beta_curr[(truth[1]+1):length(beta_curr)]
      # u transformation
      #u1 <- pu(x[,1:truth[1]], beta1)
      #u2 <- pu(x[,1:truth[2]], beta2)
      #eta1 <- u1$u
      #eta2 <- u2$u
      eta1 <- x[,1:truth[1]] %*% beta1
      eta2 <- x[,1:truth[2]] %*% beta2
      # calculate B-spline basis
      B1 <- bsplineS(eta1, breaks = quantile(eta1, knots))
      B2 <- bsplineS(eta2, breaks = quantile(eta2, knots))
      B <- cbind(z*B1, B2)
      p1 <- dim(B1)[2]
      p2 <- dim(B2)[2]
      # estimate g1 and g2
      fit.delta <- coxph(Surv(y, status) ~ 0+B)
      delta <- drop(coef(fit.delta))
      delta[is.na(delta)] <- 0
      delta1 <- delta[1 : p1]
      delta2 <- delta[(p1+1) : (p1+p2)]
      delta1[is.na(delta1)] = 0
      delta2[is.na(delta2)] = 0
      # step 2 fix g to estimate beta
      # calculate first derivative of B-spline basis
      B_deriv_1 <- bsplineS(eta1, breaks = quantile(eta1, knots), nderiv = 1)
      B_deriv_2 <- bsplineS(eta2, breaks = quantile(eta2, knots), nderiv = 1)
      # calculate input 
      #newx_1 <- z * drop((B_deriv_1*(u1$deriv))%*%delta1)*x[,1:truth[1]]
      #newx_2 <- drop((B_deriv_2*(u2$deriv))%*%delta2)*x[,1:truth[2]]
      newx_1 <- z * drop(B_deriv_1%*%delta1)*x[,1:truth[1]]
      newx_2 <- drop(B_deriv_2%*%delta2)*x[,1:truth[2]]
      newx <- cbind(newx_1, newx_2)
      # calculate offset
      off_1 <- z * (B1%*%delta1 - newx_1 %*% beta1)
      off_2 <- B2%*%delta2 - newx_2 %*% beta2
      off <- off_1 + off_2
      # estimate beta
      beta <- my_surv(newx, y, atrisk, off, lam = lam)
      #if(sum(is.na(beta)) > 1){
      #	break
      #	stop("only 1 variable in betas; decrease lambda")
      #}
      beta1 <- beta[1:truth[1]]
      beta2 <- beta[(truth[1]+1):length(beta_curr)]
      check <- c(sum(beta1!=0),sum(beta2!=0))
      if(min(check) <= 1) {
        stop("0 beta occurs; decrease lambda")
        flag <- TRUE
        if(check[1] != 0) beta[1:truth[1]] <- normalize(beta1) 
        if(check[2] !=0) beta[(truth[1]+1):length(beta_curr)] <- normalize(beta2) 
        break
      }
      beta <- c(normalize(beta1), normalize(beta2))
      conv <- (max(abs(beta - beta_curr)) < eps)
      beta_curr <- beta
    }
    geta <- z * B1 %*% delta1 + B2 %*% delta2
    g1 <- B1 %*% delta1
    g2 <- B2 %*% delta2
    loss <- -fit.delta$loglik[2]
    df1 <- sum(beta1!=0)
    df2 <- sum(beta2!=0)
    # df <- df1 + df2
    df <- df1 + df2 + 2*length(delta1)
    bic <- -2 * loss + df * log(n) * log(p)
    aic <- -2 * loss + 2 * df
    delta.var <- summary(fit.delta)$var
    out <- list(beta1 = beta[1:truth[1]], beta2 = beta[(truth[1]+1):length(beta_curr)], B1 = B1, B2 = B2, delta2 = delta2, delta1 = delta1, iter = iter, g = geta, g1 = g1, g2 = g2, loss = loss, df = df, df1 = df1, df2 = df2, bic = bic, aic = aic, x = x, y = y, z = z, status = status, knots = knots, flag = flag, p=p, conv=conv, final.x = newx, final.off = off)
    class(out) <- "cste"
    return(out)
  }
}
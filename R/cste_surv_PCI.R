
#' Calculate pointwise confidence intervals of CSTE curve for time to event outcome with right censoring.
#' 
#' This function calculates pointwise confidence intervals of CSTE curve for time to event outcome with right censoring.
#' 
#' 
#'@param fit a S3 class of cste.
#'@param alpha the pointwise confidence intervals are of \eqn{1-\alpha} confidence level.
#'
#'@return A list which includes:
#' \itemize{
#'    \item \code{or_x}: the ordered value of \eqn{X\beta_1}. 
#'    \item \code{fit_x}: the fitted value of CSTE curve corresponding to \code{or_x}.
#'    \item \code{lower_bound}: the lower bound of CSTE's pointwise confidence intervals.
#'    \item \code{upper_bound}: the upper bound of CSTE's pointwise confidence intervals.
#' }
#' @references
#' Ma Y. and Zhou X. (2017). 
#' Treatment selection in a randomized clinical trial via covariate-specific 
#' treatment effect curves, \emph{Statistical Methods in Medical Research}, 26(1), 124-141.
#' 
#' @seealso \code{\link{cste_surv}}

cste_surv_PCI <- function(fit, alpha=0.05){
  x <- fit$x
  y <- fit$y
  z <- fit$z
  status <- fit$status
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (p==1){
    u1 <- x
    u2 <- x
    h_x <- bw.nrd0(x)
    f_x <- sapply(x, function(xx) mean(dnorm((xx - x)/h_x)/h_x))
  } else {
    #u1 <- pu(x, fit$beta1)$u
    #u2 <- pu(x, fit$beta2)$u
    u1 <- x %*% fit$beta1
    u2 <- x %*% fit$beta2
    h_x <- bw.nrd0(u1)
    f_x <- sapply(u1, function(xx) mean(dnorm((xx - u1)/h_x)/h_x))
  }
  B1 <- fit$B1
  B2 <- fit$B2
  B <- cbind(z*B1, B2)
  g <- fit$g
  g1 <- fit$g1
  fit.x <- as.matrix(x %*% fit$beta1)
  fit.delta <- coxph(Surv(y,status) ~ 0+B)
  fit.delta$coefficients[is.na(fit.delta$coefficients)]=0
  fit.surv <- survfit(fit.delta, newdata = data.frame(x))
  S <- fit.surv$surv
  nevent <- fit.surv$n.event
  chaz <- basehaz(fit.delta, centered=FALSE)$hazard
  haz0 <- chaz - c(0, chaz[1:(length(chaz)-1)])
  rho <- t(S) * as.numeric(exp(g))
  if (max(rho) == 0) rho = rho + 1e-4
  a00 <- t(f_x %*% t(colMeans(rho)))
  a10 <- t(f_x %*% t(colMeans(rho*z)))
  g1_deriv_2 <- bsplineS(u1, breaks = quantile(u1, fit$knots), nderiv = 2) %*% fit$delta1
  #bias <- h_x^4 * g1_deriv_2 / 2
  bias <- 0
  Sigma <- fit.x / (colSums((a10 - a10^2/a00*haz0) * nevent))
  L <- g1 - bias - qnorm(1-alpha/2) * sqrt(Sigma/(n*h_x))
  U <- g1 - bias + qnorm(1-alpha/2) * sqrt(Sigma/(n*h_x))
  result <- cbind(fit.x, g1, L, U)
  result <- na.omit(result)
  result <- result[order(result[,1]),]
  L <- result[,3]
  U <- result[,4]
  fit.x <- result[,1]
  g1 <- result[,2]
  
  
  return(list(or_x = fit.x, fit_x = g1, 
              lower_bound = L, upper_bound = U))
}

#' Select the optimal tuning parameters in CSTE estimation for time to event outcome with right censoring.   
#'
#' select lasso penalty parameter \eqn{\lambda} for \eqn{\beta_1} and
#'\eqn{\beta_2} in CSTE estimation.
#' 
#' 
#'@param x samples of covariates which is a \eqn{n*p} matrix.
#'@param y samples of binary outcome which is a \eqn{n*1} vector.
#'@param z samples of treatment indicator which is a \eqn{n*1} vector.
#'@param status samples of censoring indicator which is a \eqn{n*1} vector. 
#' Default value is NULL, indicating no censoring.  
#'@param lam_seq a sequence for the choice of \eqn{\lambda}. 
#'@param beta_ini initial values for \eqn{(\beta_1', \beta_2')'}, default value is NULL.
#'@param nknots number of knots for the B-spline for estimating \eqn{g_1} and \eqn{g_2}.
#'@param max.iter maximum iteration for the algorithm.
#'@param eps numeric scalar \eqn{\geq} 0, the tolerance for the estimation 
#' of \eqn{\beta_1} and \eqn{\beta_2}. 
#'
#'@return A list which includes
#' \itemize{
#'    \item \code{optimal}: optimal cste within the given the sequence of \eqn{\lambda}.
#'    \item \code{bic}: BIC for the sequence of \eqn{\lambda}.
#'    \item \code{lam_seq}: the sequence of \eqn{\lambda} that is used.
#'    
#' }
#'
#' @references
#' Ma Y. and Zhou X. (2017). 
#' Treatment selection in a randomized clinical trial via covariate-specific 
#' treatment effect curves, \emph{Statistical Methods in Medical Research}, 26(1), 124-141.
#' 
#' @seealso \code{\link{cste_surv}}

select_cste_surv <- function(x, y, z, status = NULL, lam_seq, beta_ini = NULL, nknots = 1, max.iter = 2000, eps = 1e-3) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  out <- vector("list", length(lam_seq))
  beta1 <- matrix(0, p, length(lam_seq))
  beta2 <- matrix(0, p, length(lam_seq))
  for(i in 1:length(lam_seq)) {
    if(is.null(beta_ini)){
      beta_ini <- rep(normalize(rep(1, p)), 2)
    }
    if(i > 1) beta_ini <- c(out[[i-1]]$beta1, out[[i-1]]$beta2)
    out[[i]] <- cste_surv(x, y, z, status, lam = lam_seq[i], beta_ini = beta_ini, nknots = nknots, max.iter = max.iter, eps = eps)
    beta1[, i] <- out[[i]]$beta1
    beta2[, i] <- out[[i]]$beta2
    if(out[[i]]$flag | out[[i]]$df1 <= 2 | out[[i]]$df2 <= 2) {
      warnings("not all lambdas are used; decrease lambda")
      break
    }
  }
  df <- sapply(out[1:i], function(x) x$df)
  bic <- sapply(out[1:i], function(x) x$bic)
  loss <- sapply(out[1:i], function(x) x$loss)
  return(list(optimal = out[[which.min(bic)]], bic = bic, lam_seq = lam_seq[1:i], df = df, complete = out[1:i], beta1=beta1, beta2=beta2))
}

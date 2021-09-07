#' Predict the CSTE curve of new data for time to event outcome with right censoring.
#'
#' Predict the CSTE curve of new data for time to event outcome with right censoring.
#' 
#' 
#'@param obj a S3 class of cste.
#'@param newx samples of covariates which is a \eqn{m*p} matrix.
#'@param alpha (1-\eqn{\alpha}) confidence level.  
#'
#'@return A S3 class of cste, which includes
#' \itemize{
#'    \item \code{g1}: predicted \eqn{g_1(X\beta_1)}.
#'    \item \code{lower_bound}: the lower bound of CSTE's pointwise confidence intervals.
#'    \item \code{upper_bound}: the upper bound of CSTE's pointwise confidence intervals.
#' }
#'
#' @references
#' Ma Y. and Zhou X. (2017). 
#' Treatment selection in a randomized clinical trial via covariate-specific 
#' treatment effect curves, \emph{Statistical Methods in Medical Research}, 26(1), 124-141.
#' 
#' @seealso \code{\link{cste_surv}}

predict_cste_surv <- function(obj, newx, alpha = 0.05) {
  #type <- match.arg(type)
  if(missing(newx)) {
    out <- as.numeric(obj$B1 %*% obj$delta1)
    confidence_interval <- cste_surv_PCI(obj, alpha = 0.05)
    L <- as.numeric(confidence_interval$lower_bound)
    U <- as.numeric(confidence_interval$upper_bound)
  } else {
    # if(type == "x") 
    eta <- as.numeric(newx %*% obj$beta1)
    # else eta <- as.numeric(newx)
    eta.fit <- as.numeric(obj$x %*% obj$beta1)
    matched <- sapply(eta, function(xx) which.min(abs(xx-eta.fit)))
    confidence_interval <- cste_surv_PCI(obj, alpha = 0.05)
    out <- as.numeric(confidence_interval$fit_x[matched,])
    L <- as.numeric(confidence_interval$lower_bound[matched])
    U <- as.numeric(confidence_interval$upper_bound[matched])
  }
  return(list(g1 = out, lower_bound = L, upper_bound = U))
}
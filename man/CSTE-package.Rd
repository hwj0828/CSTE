\name{CSTE-package}
\alias{CSTE-package}
\alias{CSTE}
\docType{package}
\title{
  Covariate specific treatment effect (CSTE) curve. 
}
\description{
  Estimation of optimal individualized treatment rules using a covariate specific treatment effect (CSTE) curve.  
  Suppose we have a binary treatment \eqn{Z}, a \eqn{p}-dimensional covariate \eqn{X},
  which may be high-dimensional. 
 
 If the outcome \eqn{Y} is a binary variable, we consider the CSTE curve, defined by     \deqn{ CSTE(x)  = logit( E(Y(1)\mid X=x)) - logit(E(Y(0)\mid X=x)), }
where \eqn{logit(u) = \log(u) -\log(1-u)}, \eqn{Y(1)} and \eqn{Y(0)} denote the potential outcomes if the active treatment and the control treatment are received, respectively. In this case, the CSTE curve represents the difference of the logarithm of odds ratio between treated group and control group, which is a common causal quantity in clinical studies for binary outcome.  Moreover,  
\eqn{Y = Z Y(1) + (1-Z) Y(0)}. Under the unconfoundedness assumption such that
 \eqn{(Y(0), Y(1)) \perp Z \mid  X}, the CSTE curve can be re-expressed as 
    \deqn{ CSTE(x) = logit(E(Y \mid  X=x, Z=1)) - logit( E(Y \mid  X=x, Z=1)).} 
Denoting \eqn{\mu(X, Z) = E(Y| X, Z)}, we model the logarithm of odds ration as
 \deqn{logit(\mu(X, Z)) = g_1(X\beta_1)Z + g_2(X\beta_2),} 
 Then we have
  \deqn{ CSTE(x) = g_1(x \beta_1).} 
  
  If the outcome \eqn{T} is a failure time with right censoring, we let 
\eqn{T(z)} be the failure time under the corresponding treatment arm \eqn{Z=z} 
  (\eqn{z=0,1}), and define the potential conditional hazard rate function as
\deqn{
\lambda^{(z)}(t\mid x) = \lim_{\Delta t\to 0}\frac{P(t<T(z)\leq t+\Delta t\mid T(z)>t, X=x)}{\Delta t}.}
In most medical research with a survival outcome, the treatment effect is usually represented by
the logarithm of a hazard ratio. Hence, the CSTE curve at a fixed time \eqn{t} is defined as 
  \deqn{
CSTE(x,t) = \log\left\{\frac{\lambda^{(1)}(t\mid x)}{\lambda^{(0)}(t\mid x)}\right\}.
   }
We focus on a special case that \eqn{CSTE(x,t) = CSTE(x)} in a 
randomized clinical trial. This assumption holds under many well-known models, such as Cox’s model and varying-coefficient Cox’s model. 
Consequently, we assume the following varying coefficient proportional hazard regression model for event time
  \deqn{\lambda(t| X, Z) = \lambda_0(t) \exp(g_1(X\beta_1)Z + g_2(X\beta_2)),}
which implies that 
  \deqn{   CSTE(x) = g_1(x \beta_1) }  

  This package facilitate estimating, making inference and predicting \eqn{CSTE(x)}
for binary outcome and time to event outcome with right censoring.  
}

\details{
  The R package CSTE - version 1.0 can be used for three main tasks:
  \itemize{
   \item \code{\link{cste_bin}}: estimate the CSTE curve for binary outcome.
   \item \code{\link{cste_surv}}: estimate the CSTE curve for time to event outcome 
   with right censoring.
   \item \code{\link{predict_cste_bin}}: predict the CSTE curve of new subject with binary outcome.
   \item \code{\link{predict_cste_surv}}: predict the CSTE curve of new subject for time to event outcome with right censoring.
   \item \code{\link{cste_bin_SCB}}: compute the simultaneous confidence bands (SCB) of CSTE curve for binary outcome.
   \item \code{\link{cste_surv_PCI}}: compute the point confidence intervals (PCI) of CSTE curve for time to event outcome with right censoring. 
}
}
\author{
Peng Wu, Wenjie Hu, Yuhao Deng, Xiaohua Zhou.

Maintainer: Wenjie Hu <huwenjie@pku.edu.cn>
}
\references{
Guo W., Zhou X. and Ma S. (2021). Estimation of Optimal Individualized Treatment Rules
Using a Covariate-Specific Treatment Effect Curve with High-dimensional Covariates,
 \emph{Journal of the American Statistical Association}, 116(533), 309-321.
 
 Ma Y. and Zhou X. (2017). Treatment selection in a randomized clinical trial via covariate-specific treatment effect curves, \emph{Statistical Methods in Medical Research}, 26(1), 124-141.
}
\keyword{ package }

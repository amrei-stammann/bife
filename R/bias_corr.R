#' @title
#' Asymptotic bias-correction for binary choice Models with fixed effects
#' @description
#' \code{\link{bias_corr}} is a post-estimation routine that can be used to substantially reduce the 
#' incidental parameter bias problem (Neyman and Scott (1948)) present in non-linear fixed effects 
#' models (see Fernandez-Val and Weidner (2018) for an overview). The command applies the analytical 
#' bias-correction derived by Fernandez-Val (2009) to obtain bias-corrected estimates of 
#' the structural parameters.
#' 
#' \strong{Remark:} Fernandez-Val (2009) further refined the bias-correction of 
#' Hahn and Newey (2004). The correction is now also applicable to dynamic models.
#' @param
#' object an object of class \code{"bife"}.
#' @param
#' L unsigned integer indicating a bandwidth for the estimation of spectral densities proposed by 
#' Hahn and Kuersteiner (2011). Default is zero, which should be used if all regressors are 
#' assumed to be strictly exogenous. In the presence of weakly exogenous or predetermined 
#' regressors, Fernandez-Val and Weidner (2018) suggest to choose a bandwidth not higher than 
#' four.
#' @return
#' The function \code{\link{bias_corr}} returns a named list of class \code{"bife"}.
#' @references
#' Fernandez-Val, I. (2009). "Fixed effects estimation of structural parameters and marginal 
#' effects in panel probit models". Journal of Econometrics 150(1), 71-85.
#' @references
#' Fernandez-Val, I. and Weidner, M. (2018). "Fixed effects estimation of large-t panel data 
#' models". Annual Review of Economics, 10, 109-138.
#' @references
#' Hahn, J. and Kuersteiner, G. (2011). "Bias reduction for dynamic nonlinear panel models with 
#' fixed effects". Econometric Theory, 27(6), 1152-1191.
#' @references 
#' Hahn, J. and Newey, W. (2004). "Jackknife and analytical bias reduction for nonlinear panel 
#' models". Econometrica 72(4), 1295-1319.
#' @references
#' Neyman, J. and Scott, E. L. (1948). "Consistent estimates based on partially consistent 
#' observations". Econometrica, 16(1), 1-32.
#' @references
#' Stammann, A., Heiss, F., and and McFadden, D. (2016). "Estimating Fixed Effects Logit Models with 
#' Large Panel Data". Working paper.
#' @seealso
#' \code{\link{bife}}
#' @examples 
#' \donttest{
#' # Load 'psid' dataset
#' library(bife)
#' dataset <- psid
#' 
#' # Fit a static logit model
#' mod <- bife(LFP ~ I(AGE^2) + log(INCH) + KID1 + KID2 + KID3 + factor(TIME) | ID, dataset)
#' summary(mod)
#' 
#' # Apply analytical bias-correction
#' mod_bc <- bias_corr(mod)
#' summary(mod_bc)
#' }
#' @export
bias_corr <- function(object, L = 0L) {
  # Validity of input argument (object)
  if (!inherits(object, "bife")) {
    stop("'bias_corr' called on a non-'bife' object.", call. = FALSE)
  }
  
  # Information in case of predetermined variables
  if (L > 0L) {
    message("Ensure that your data set is sorted by time!")
  }
  
  # Get names of the fixed effects variables and extract individual number of time periods
  idvar <- attr(terms(object[["formula"]], rhs = 2L), "term.labels")
  Ti <- object[["data"]][, .N, by = eval(idvar)][[2L]]
  
  # Extract model response and regressor matrix
  y <- object[["data"]][[1L]]
  X <- model.matrix(object[["formula"]], object[["data"]], rhs = 1L)[, - 1L, drop = FALSE]
  id <- as.integer(object[["data"]][[idvar]])
  attr(X, "dimnames") <- NULL
  
  # Compute required derivatives
  beta <- object[["coefficients"]]
  family <- object[["family"]]
  eta <- as.vector(X %*% beta + object[["alpha"]][id])
  mu <- family[["linkinv"]](eta)
  mu_eta <- family[["mu.eta"]](eta)
  if (family[["link"]] == "logit") {
    v <- y - mu
    w <- mu_eta
    z <- w * (1.0 - 2.0 * mu)
  } else {
    w <- mu_eta / family[["variance"]](mu)
    v <- w * (y - mu)
    w <- w * mu_eta
    z <- - eta * w
  }
  
  # Bias correction of Fernandez-Val (2009)
  MX <- center_variables(X * sqrt(w), sqrt(w), Ti) / sqrt(w)
  B <- as.vector(group_sums_bias(MX * z, w, Ti)) / 2.0
  if (L > 0L) {
    B <- B + as.vector(group_sums_spectral(MX * w, v, w, L, Ti))
  }
  beta <- beta - solve(object[["Hessian"]], - B)
  rm(mu, v, MX)
  
  # Adjust \alpha using an offset algorithm and update \eta
  eta <- as.vector(X %*% beta)
  alpha <- bife_offset(y, eta, id, Ti, family, object[["control"]])
  
  # Recompute weights
  eta <- eta + alpha[id]
  mu_eta <- family[["mu.eta"]](eta)
  if (family[["link"]] == "logit") {
    w <- mu_eta
  } else {
    w <- mu_eta^2 / family[["variance"]](family[["linkinv"]](eta))
  }
  
  # Recompute Hessian
  H <- crossprod(center_variables(X * sqrt(w), sqrt(w), Ti))
  
  # Recompute deviance
  mu <- family[["linkinv"]](eta)
  dev <- sum(family[["dev.resids"]](y, mu, rep(1.0, length(y))))
  
  # Add names and modify result list
  dimnames(H) <- list(names(beta), names(beta))
  names(alpha) <- names(object[["alpha"]])
  object[["coefficients"]] <- beta
  object[["alpha"]] <- alpha
  object[["deviance"]] <- dev
  object[["Hessian"]] <- H
  object[["bandwidth"]] <- L
  
  # Return updated result list
  object
}
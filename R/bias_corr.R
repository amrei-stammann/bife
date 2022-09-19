#' @title
#' Asymptotic bias correction for binary choice Models with fixed effects
#' @description
#' \code{\link{bias_corr}} is a post-estimation routine that can be used to substantially reduce the 
#' incidental parameter bias problem (Neyman and Scott (1948)) present in non-linear fixed effects 
#' models (see Fernández-Val and Weidner (2018) for an overview). The command applies the analytical 
#' bias correction derived by Fernández-Val (2009) to obtain bias-corrected estimates of 
#' the structural parameters.
#' 
#' \strong{Remark:} Fernández-Val (2009) further refined the bias correction of 
#' Hahn and Newey (2004). The correction is now also applicable to models with weakly exogenous regressors.
#' @param
#' object an object of class \code{"bife"}.
#' @param
#' L unsigned integer indicating a bandwidth for the estimation of spectral densities proposed by 
#' Hahn and Kuersteiner (2011). Default is zero, which should be used if all regressors are 
#' assumed to be strictly exogenous. In the presence of weakly exogenous or predetermined 
#' regressors, Fernández-Val and Weidner (2018) suggest to choose a bandwidth not higher than 
#' four.
#' @return
#' The function \code{\link{bias_corr}} returns a named list of class \code{"bife"}.
#' @references
#' Fernández-Val, I. (2009). "Fixed effects estimation of structural parameters and marginal 
#' effects in panel probit models". Journal of Econometrics 150(1), 71-85.
#' @references
#' Fernández-Val, I. and M. Weidner (2018). "Fixed effects estimation of large-t panel data 
#' models". Annual Review of Economics, 10, 109-138.
#' @references
#' Hahn, J. and G. Kuersteiner (2011). "Bias reduction for dynamic nonlinear panel models with 
#' fixed effects". Econometric Theory, 27(6), 1152-1191.
#' @references 
#' Hahn, J. and W. Newey (2004). "Jackknife and analytical bias reduction for nonlinear panel 
#' models". Econometrica 72(4), 1295-1319.
#' @references
#' Neyman, J. and E. L. Scott (1948). "Consistent estimates based on partially consistent 
#' observations". Econometrica, 16(1), 1-32.
#' @references
#' Stammann, A., F. Heiss, and D. McFadden (2016). "Estimating Fixed Effects Logit Models with 
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
#' # Apply analytical bias correction
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
  beta_uncorr <- object[["coefficients"]]
  alpha_uncorr <- object[["alpha"]]
  family <- object[["family"]]
  eta <- as.vector(X %*% beta_uncorr) + alpha_uncorr[id]
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
  
  # Bias correction of Fernández-Val (2009)
  nt <- object[["nobs"]][["nobs"]]
  MX <- center_variables(X * sqrt(w), sqrt(w), Ti) / sqrt(w)
  B <- as.vector(group_sums_bias(MX * z, w, Ti)) / 2.0 / nt
  if (L > 0L) {
    B <- B + as.vector(group_sums_spectral(MX * w, v, w, L, Ti)) / nt
  }
  bias_term <- solve(object[["Hessian"]] / nt, - B)
  beta <- beta_uncorr - bias_term
  rm(mu, v, MX)
  
  # Adjust \alpha using an offset algorithm and update the linear index
  xb <- as.vector(X %*% beta)
  alpha <- bife_offset(y, xb, id, Ti, family, object[["control"]])
  
  # Recompute Hessian
  eta <- xb + alpha[id]
  mu_eta <- family[["mu.eta"]](eta)
  if (family[["link"]] == "logit") {
    w <- mu_eta
  } else {
    w <- mu_eta^2 / family[["variance"]](family[["linkinv"]](eta))
  }
  H <- crossprod(center_variables(X * sqrt(w), sqrt(w), Ti))
  
  # Add names and modify result list
  names(beta) <- names(beta_uncorr)
  names(alpha) <- names(alpha_uncorr)
  dimnames(H) <- list(names(beta_uncorr), names(beta_uncorr))
  object[["coefficients"]] <- beta
  object[["coefficients_uncorr"]] <- beta_uncorr
  object[["alpha"]] <- alpha
  object[["alpha_uncorr"]] <- alpha_uncorr
  object[["bias_term"]] <- bias_term
  object[["bandwidth"]] <- L
  object[["Hessian"]] <- H
  
  # Return updated result list
  object
}
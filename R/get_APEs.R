#' @title
#' Compute average partial effects for binary choice models with fixed effects
#' @description
#' \code{\link{get_APEs}} is a post-estimation routine that can be used to estimate average partial 
#' effects with respect to all covariates in the model and the corresponding covariance matrix. The 
#' estimation of the covariance is based on a linear approximation (delta method). Note that 
#' the command automatically determines which of the regressors are continuous or binary.
#' 
#' \strong{Remark:} The routine currently does not allow to compute average partial effects based 
#' on functional forms like interactions and polynomials.
#' 
#' \strong{Note:} \code{\link{apeff_bife}} is deprecated and will be removed soon.
#' @param
#' object an object of class \code{"bife"}.
#' @param
#' n_pop unsigned integer indicating a finite population correction for the estimation of the 
#' covariance matrix of the average partial effects proposed by 
#' Cruz-Gonzalez, Fernández-Val, and Weidner (2017). The correction factor is computed as follows: 
#' \eqn{(n^{\ast} - n) / (n^{\ast} - 1)}{(n.pop - n) / (n.pop - 1)}, 
#' where \eqn{n^{\ast}}{n.pop} and \eqn{n}{n} are the size of the entire 
#' population and the full sample size. Default is \code{NULL}, which refers to a factor of one and 
#' is equal to an infinitely large population.
#' @param
#' sampling_fe a string equal to \code{"independence"} or \code{"unrestricted"} which imposes 
#' sampling assumptions about the unobserved effects. \code{"independence"} imposes that all 
#' unobserved effects are mutually independent sequences. \code{"unrestricted"} does not impose any
#' sampling assumptions. Note that this option only affects the estimation of the covariance. 
#' Default is \code{"independence"}.
#' @param
#' weak_exo logical indicating if some of the regressors are assumed to be weakly exogenous (e.g. 
#' predetermined). If object is returned by \code{\link{bias_corr}}, the option will be 
#' automatically set to \code{TRUE} if the choosen bandwidth parameter is larger than zero. 
#' Note that this option only affects the estimation of the covariance matrix. 
#' Default is \code{FALSE}, which assumes that all regressors are strictly exogenous.
#' @param
#' ... arguments passed to the deprecated function \code{\link{apeff_bife}}.
#' @return
#' The function \code{\link{get_APEs}} returns a named list of class \code{"bifeAPEs"}.
#' @references
#' Cruz-Gonzalez, M., I. Fernández-Val, and M. Weidner. (2017). "Bias corrections for probit and 
#' logit models with two-way fixed effects". The Stata Journal, 17(3), 517-545.
#' @references
#' Fernández-Val, I. (2009). "Fixed effects estimation of structural parameters and marginal 
#' effects in panel probit models". Journal of Econometrics 150(1), 71-85.
#' @references
#' Fernández-Val, I. and M. Weidner (2018). "Fixed effects estimation of large-t panel data models". 
#' Annual Review of Economics, 10, 109-138.
#' @references
#' Neyman, J. and E. L. Scott (1948). "Consistent estimates based on partially consistent 
#' observations". Econometrica, 16(1), 1-32.
#' @references
#' Stammann, A., F. Heiss, and D. McFadden (2016). "Estimating Fixed Effects Logit Models with 
#' Large Panel Data". Working paper.
#' @seealso
#' \code{\link{bias_corr}}, \code{\link{bife}}
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
#' # Compute average partial effects
#' mod_ape <- get_APEs(mod)
#' summary(mod_ape)
#' 
#' # Apply analytical bias correction
#' mod_bc <- bias_corr(mod)
#' summary(mod_bc)
#' 
#' # Compute bias-corrected average partial effects
#' mod_ape_bc <- get_APEs(mod_bc)
#' summary(mod_ape_bc)
#' }
#' @export
get_APEs <- function(object,
                     n_pop       = NULL,
                     sampling_fe = c("independence", "unrestricted"),
                     weak_exo    = FALSE) {
  # Validity of input argument (object)
  if (!inherits(object, "bife")) {
    stop("'get_APEs' called on a non-'bife' object.", call. = FALSE)
  }
  
  # Check validity of 'sampling_fe'
  sampling_fe <- match.arg(sampling_fe)
  
  # Extract prior information if available
  biascorr <- !is.null(object[["bandwidth"]])
  if (biascorr) {
    L <- object[["bandwidth"]]
    if (L > 0L) {
      weak_exo <- TRUE
    } else {
      weak_exo <- FALSE
    }
  }
  
  # Get names of the fixed effects variables and extract individual number of time periods
  idvar <- attr(terms(object[["formula"]], rhs = 2L), "term.labels")
  Ti <- object[["data"]][, .N, by = eval(idvar)][[2L]]
  
  # Extract model information
  beta <- object[["coefficients"]]
  family <- object[["family"]]
  nt_full <- object[["nobs"]][["nobs_full"]]
  nt <- object[["nobs"]][["nobs"]]
  p <- length(beta)
  
  # Check validity of 'n_pop'
  if (!is.null(n_pop)) {
    n_pop <- as.integer(n_pop)
    if (n_pop < nt_full) {
      warning(paste("Size of the entire population is lower than the full sample size.",
                    "Correction factor set to zero."), call. = FALSE)
      adj <- 0.0
    } else {
      adj <- (n_pop - nt_full) / (n_pop - 1L)
    }
  } else {
    adj <- 1.0
  }
  
  # Extract model response and regressor matrix
  y <- object[["data"]][[1L]]
  X <- model.matrix(object[["formula"]], object[["data"]], rhs = 1L)[, - 1L, drop = FALSE]
  id <- as.integer(object[["data"]][[idvar]])
  attr(X, "dimnames") <- NULL
  
  # Determine which of the regressors are binary
  binary <- apply(X, 2L, function(x) all(x %in% c(0L, 1L)))
  
  # Compute required derivatives
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
  rm(y)
  
  # Compute average partial effects and Jacobian
  PX <- X - center_variables(X * sqrt(w), sqrt(w), Ti) / sqrt(w)
  Delta <- matrix(NA_real_, nt, p)
  Delta1 <- matrix(NA_real_, nt, p)
  J <- matrix(NA_real_, p, p)
  Delta[, !binary] <- mu_eta
  Delta1[, !binary] <- partial_mu_eta(eta, family, 2L)
  for (i in seq.int(p)) {
    if (binary[[i]]) {
      eta0 <- eta - X[, i] * beta[[i]]
      f1 <- family[["mu.eta"]](eta0 + beta[[i]])
      Delta[, i] <- family[["linkinv"]](eta0 + beta[[i]]) - family[["linkinv"]](eta0)
      Delta1[, i] <- f1 - family[["mu.eta"]](eta0)
      J[i, ] <- - colSums(PX * Delta1[, i])
      J[i, i] <- sum(f1) + J[i, i]
      J[i, - i] <- colSums(X[, - i, drop = FALSE] * Delta1[, i]) + J[i, - i]
      rm(eta0, f1)
    } else {
      Delta[, i] <- beta[[i]] * Delta[, i]
      Delta1[, i] <- beta[[i]] * Delta1[, i]
      J[i, ] <- colSums((X - PX) * Delta1[, i])
      J[i, i] <- sum(mu_eta) + J[i, i]
    }
  }
  delta <- colSums(Delta) / nt_full
  Delta <- t(t(Delta) - delta)
  rm(mu_eta)
  
  # Compute projection of \Psi
  Psi <- - Delta1 / w
  MPsi <- center_variables(Psi * sqrt(w), sqrt(w), Ti) / sqrt(w)
  PPsi <- Psi - MPsi
  rm(Delta1, Psi)
  
  # Compute analytical bias correction of average partial effects
  if (biascorr) {
    # Compute second-order partial derivatives
    Delta2 <- matrix(NA_real_, nt, p)
    Delta2[, !binary] <- partial_mu_eta(eta, family, 3L)
    for (i in seq.int(p)) {
      if (binary[[i]]) {
        eta0 <- eta - X[, i] * beta[[i]]
        Delta2[, i] <- partial_mu_eta(eta0 + beta[[i]], family, 2L) - 
          partial_mu_eta(eta0, family, 2L)
        rm(eta0)
      } else {
        Delta2[, i] <- beta[[i]] * Delta2[, i]
      }
    }
    
    # Compute corrected average partial effect using the bias correction of Fernández-Val (2009)
    B <- as.vector(group_sums_bias(PPsi * z + Delta2, w, Ti)) / 2.0
    rm(Delta2)
    if (L > 0L) {
      B <- B - as.vector(group_sums_spectral(MPsi * w, v, w, L, Ti))
    }
    delta <- delta - B / nt_full
  }
  rm(eta, mu, MPsi)
  
  # Compute standard errors
  J <- J %*% solve(object[["Hessian"]])
  Gamma <- tcrossprod((X - PX) * v, J) - PPsi * v
  V <- crossprod(Gamma)
  if (adj > 0.0) {
    # Simplify covariance if sampling assumptions are imposed
    if (sampling_fe == "independence") {
      V <- V + adj * group_sums_var(Delta, Ti)
    } else {
      V <- V + adj * tcrossprod(colSums(Delta))
    }
    
    # Add covariance in case of weak exogeneity
    if (weak_exo) {
      V <- V + adj * 2.0 * group_sums_cov(Delta, Gamma, Ti)
    }
  }
  V <- V / nt_full^2
  
  # Add names
  names(delta) <- names(beta)
  dimnames(V) <- list(names(beta), names(beta))
  
  # Return result list
  structure(list(delta       = delta,
                 vcov        = V,
                 sampling_fe = sampling_fe,
                 weak_exo    = weak_exo),
            class = "bifeAPEs")
}


### Deprecated functions

#' @rdname get_APEs
#' @aliases get_APEs
#' @export
apeff_bife <- function(...) {
  .Deprecated("get_APEs")
}
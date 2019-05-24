#' @title
#' Efficiently fit binary choice models with fixed effects
#' @description
#' \code{\link{bife}} can be used to fit fixed effects binary choice models (logit and probit) 
#' based on an unconditional maximum likelihood approach. It is tailored for the fast estimation of 
#' binary choice models with potentially many individual fixed effects. The routine is based on a 
#' special pseudo demeaning algorithm derived by Stammann, Heiss, and McFadden (2016). The 
#' estimates obtained are identical to the ones of \code{\link[stats]{glm}}, but the computation 
#' time of \code{\link{bife}} is much lower.
#' 
#' \strong{Remark:} The term fixed effect is used in econometrician's sense of having a full set of 
#' individual specific intercepts. All other parameters in the model are referred to as 
#' structural parameters.
#' @param
#' formula an object of class \code{"formula"} (or one that can be coerced to that class): 
#' a symbolic description of the model to be fitted. \code{formula} must be of type 
#' \eqn{y ~ x | id} where the \code{id} refers to an individual identifier (fixed effect category).
#' @param
#' data an object of class \code{"data.frame"} containing the variables in the model.
#' @param 
#' model the description of the error distribution and link function to be used in the model. 
#' For \code{\link{bife}} this has to be a character string naming the model function. 
#' Default is \code{"logit"}.
#' @param 
#' beta_start an optional vector of starting values used for the structural parameters in the 
#' optimization algorithm. Default is zero for all structural parameters.
#' @param
#' control a named list of parameters for controlling the fitting process. See 
#' \code{\link{bife_control}} for details.
#' @param
#' bias_corr deprecated; see \code{\link{bias_corr}}.
#' @param 
#' tol_demeaning,iter_demeaning,tol_offset,iter_offset deprecated; see \code{\link{bife_control}}.
#' @details
#' \code{\link{bife}} drops all observations of cross-sectional units (individuals) with 
#' non-varying response. This can de done because these observations do not contribute to the 
#' identification of the structural parameters (perfect classification).
#' 
#' If \code{\link{bife}} does not converge this is usually a sign of linear dependence between 
#' one or more regressors and the fixed effects. In this case, you should carefully inspect 
#' your model specification.
#' @return
#' The function \code{\link{bife}} returns a named list of class \code{"bife"}.
#' @references
#' Stammann, A., Heiss, F., and and McFadden, D. (2016). "Estimating Fixed Effects Logit Models with 
#' Large Panel Data". Working paper.
#' @examples 
#' \donttest{
#' # Load 'psid' dataset
#' library(bife)
#' dataset <- psid
#' 
#' # Fit a static logit model
#' mod <- bife(LFP ~ I(AGE^2) + log(INCH) + KID1 + KID2 + KID3 + factor(TIME) | ID, dataset)
#' summary(mod)
#' }
#' @importFrom data.table setDT first setattr setkeyv .N .SD :=
#' @importFrom Formula Formula
#' @importFrom stats binomial coef model.frame model.matrix pnorm printCoefmat terms vcov
#' @importFrom Rcpp evalCpp
#' @useDynLib bife, .registration = TRUE 
#' @export
bife <- function(formula, data = list(), model = c("logit", "probit"),
                 beta_start = NULL, control = list(),
                 bias_corr = NULL, tol_demeaning = NULL, iter_demeaning = NULL,
                 tol_offset = NULL, iter_offset = NULL) {
  # Notification that bias-corrections are now a post-estimation routine
  if (!is.null(bias_corr)) {
    warning("Bias-correction is transfered to the post-estimation routine 'bias_corr'.",
            call. = FALSE)
  }
  
  # Deprecated function arguments - Will be removed soon!
  if (!is.null(tol_demeaning) || !is.null(iter_demeaning) ||
      !is.null(tol_offset) || !is.null(iter_offset)) {
    warning(paste0("Some function arguments are deprecated; ",
                   "please see 'bife_control' and use 'control' instead."), call. = FALSE)
  }
  
  # Validity of input argument (model)
  model <- match.arg(model)
  family <- binomial(model)
  
  # Validity of input argument (control)
  if (!inherits(control, "list")) {
    stop("'control' has to be of class list.", call. = FALSE)
  }
  
  # Extract control list
  control <- do.call(bife_control, control)
  
  # Update formula and do further validity check
  formula <- Formula(formula)
  if (length(formula)[[2L]] != 2L || length(formula)[[1L]] > 1L) {
    stop("'formula' uncorrectly specified.", call. = FALSE)
  }
  
  # Generate model.frame
  data <- suppressWarnings(model.frame(formula, data))
  setDT(data)
  lhs <- names(data)[[1L]]
  nobs_full <- nrow(data)
  nobs_na <- length(attr(data, "na.action"))
  
  # Ensure that model response is in line with the choosen model
  if (data[, is.numeric(get(lhs))]) {
    # Check if 'y' is in [0, 1]
    if (data[, any(get(lhs) < 0.0 | get(lhs) > 1.0)]) {
      stop("Model response has to be within the unit interval.", call. = FALSE)
    }
  } else {
    # Check if 'y' is factor and transform otherwise
    data[, (1L) := check_factor(get(lhs))]
    
    # Check if the number of levels equals two
    if (data[, length(levels(get(lhs)))] != 2L) {
      stop("Model response has to be binary.", call. = FALSE)
    }
    
    # Ensure 'y' is 0-1 encoded
    data[, (1L) := as.integer(get(lhs)) - 1L]
  }
  
  # Get names of the fixed effects variables and sort data
  idvar <- attr(terms(formula, rhs = 2L), "term.labels")
  setkeyv(data, idvar)
  
  # Drop observations that do not contribute to the loglikelihood
  trms <- attr(data, "terms") # Store terms; required for model matrix
  tmpvar <- temp_var(data)
  data[, (tmpvar) := mean(get(lhs)), by = eval(idvar)]
  data <- data[get(tmpvar) > 0.0 & get(tmpvar) < 1.0]
  setattr(data, "terms", trms)
  rm(trms)
  
  # Transform fixed effects indicator to factor and generate auxiliary variables
  data[, (idvar) := lapply(.SD, check_factor), .SDcols = idvar]
  nms_fe <- data[, levels(get(idvar))]
  Ti <- data[, .N, by = eval(idvar)][[2L]]

  # Determine number of dropped observations
  nobs <- c(nobs_full = nobs_full,
            nobs_na   = nobs_na,
            nobs_pc   = nobs_full - nrow(data),
            nobs      = nrow(data))
  
  # Extract model response and regressor matrix
  y <- data[[1L]]
  X <- model.matrix(formula, data, rhs = 1L)[, - 1L, drop = FALSE]
  id <- as.integer(data[[idvar]])
  nms_sp <- attr(X, "dimnames")[[2L]] # Saves memory
  attr(X, "dimnames") <- NULL

  # Check for linear dependence in 'X'
  p <- ncol(X)
  if (qr(X)[["rank"]] < p) {
    stop("Linear dependent terms detected!", call. = FALSE)
  }
  
  # Compute or check validity of starting guesses
  if (is.null(beta_start)) {
    beta <- numeric(p)
  } else {
    # Check validity of starting guesses
    if (length(beta_start) != p) {
      stop("'beta_start' must be of same dimension as the number of structural parameters.",
           call. = FALSE)
    }
    beta <- beta_start
  }
  rm(beta_start)
  
  # Compute starting guesses
  alpha <- family[["linkfun"]]((data[, first(get(tmpvar)), by = eval(idvar)][[2L]] + 0.5) / 2.0)
  data[, (tmpvar) := NULL]
  
  # Fit bife
  fit <- bife_fit(beta, alpha, y, X, id, Ti, family, control)
  
  # Add names to \beta, \alpha, and Hessian
  names(fit[["coefficients"]]) <- nms_sp
  names(fit[["alpha"]]) <- nms_fe
  dimnames(fit[["Hessian"]]) <- list(nms_sp, nms_sp)
  
  # Return result list
  structure(c(fit, list(nobs     = nobs,
                        formula  = formula,
                        data     = data,
                        family   = family,
                        control  = control)), class = "bife")
}


#' @title
#' Set \code{bife} Control Parameters
#' @description
#' Set and change parameters used for fitting \code{\link{bife}}.
#' @param
#' dev_tol tolerance level for the first stopping condition of the maximization routine. The 
#' stopping condition is based on the relative change of the deviance in iteration \eqn{r}
#' and can be expressed as follows: \eqn{(dev_{r - 1} - dev_{r}) / (0.1 + dev_{r}) < 
#' tol}{\Delta dev / (0.1 + dev) < tol}. Default is \code{1.0e-08}.
#' @param
#' rho_tol tolerance level for the stephalving in the maximization routine. Stephalving only takes
#' place if the deviance in iteration \eqn{r} is larger than the one of the previous iteration. If 
#' this is the case, 
#' \eqn{||\boldsymbol{\beta}_{r} - \boldsymbol{\beta}_{r - 1}||_{2}}{||\Delta \beta||} is 
#' halfed until the deviance is less or numerically equal compared to the deviance of the previous
#' iteration. Stephalving fails if the the following condition holds: \eqn{\rho < tol}{\rho < tol}, 
#' where \eqn{\rho}{\rho} is the stepcorrection factor. If stephalving fails the maximization
#' routine is canceled. Default is \code{1.0e-04}.
#' @param
#' conv_tol tolerance level that accounts for rounding errors inside the stephalving routine when
#' comparing the deviance with the one of the previous iteration. Default is \code{1.0e-06}.
#' @param
#' iter_max unsigned integer indicating the maximum number of iterations in the maximization
#' routine. Default is \code{100L}.
#' @param
#' trace logical indicating if output should be produced in each iteration. Default is \code{FALSE}.
#' @return
#' The function \code{\link{bife_control}} returns a named list of control 
#' parameters.
#' @seealso
#' \code{\link{bife}}
#' @export
bife_control <- function(dev_tol        = 1.0e-08,
                         rho_tol        = 1.0e-04,
                         conv_tol       = 1.0e-06,
                         iter_max       = 100L,
                         trace          = FALSE) {
  # Check validity of tolerance parameters
  if (dev_tol <= 0.0 || rho_tol <= 0.0 || conv_tol <= 0.0) {
    stop("All tolerance paramerters should be greater than zero.", call. = FALSE)
  }
  
  # Check validity of 'iter_max'
  iter_max <- as.integer(iter_max)
  if (iter_max < 1L) {
    stop("Maximum number of iterations should be at least one.", call. = FALSE)
  }
  
  # Return list with control parameters
  list(dev_tol  = dev_tol,
       rho_tol  = rho_tol,
       conv_tol = conv_tol,
       iter_max = iter_max,
       trace    = as.logical(trace))
}


### Internal functions (not exported)


# Fitting function
bife_fit <- function(beta, alpha, y, X, id, Ti, family, control) {
  # Compute initial quantities for the maximization routine
  eta <- as.vector(X %*% beta + alpha[id])
  mu <- family[["linkinv"]](eta)
  wt <- rep(1.0, length(y))
  dev <- sum(family[["dev.resids"]](y, mu, wt))
  
  # Start maximization of the log-likelihood
  conv <- FALSE
  for (iter in seq.int(control[["iter_max"]])) {
    # Compute weights and dependent variable
    mu_eta <- family[["mu.eta"]](eta)
    w_tilde <- sqrt(mu_eta^2 / family[["variance"]](mu))
    nu <- (y - mu) / mu_eta
    
    # Center regressor matrix
    MX <- center_variables(X * w_tilde, w_tilde, Ti)
    
    # Compute update steps
    beta_upd <- qr.solve(MX, nu * w_tilde)
    alpha_upd <- as.vector(update_alpha(as.vector(nu - X %*% beta_upd) * w_tilde, w_tilde, Ti))
    
    # Step-halving based on residual deviance as common in glm's
    dev_old <- dev
    rho <- 1.0
    repeat {
      # Compute residual deviance
      eta <- as.vector(X %*% (beta + rho * beta_upd) + (alpha + rho * alpha_upd)[id])
      mu <- family[["linkinv"]](eta)
      dev <- sum(family[["dev.resids"]](y, mu, wt))
      
      # Check if deviance is not increasing
      if (is.finite(dev) && dev <= dev_old + control[["conv_tol"]] * dev) {
        # Update \beta and \alpha
        beta <- beta + rho * beta_upd
        alpha <- alpha + rho * alpha_upd
        break
      }
      
      # Update \rho
      rho <- rho / 2.0
      
      # If \rho is to small throw error
      if (rho < control[["rho_tol"]]) {
        stop("Failure in step-halving.", call. = FALSE)
      }
    }
    
    # Progress information
    if (control[["trace"]]) {
      cat("Deviance=", dev, "- Iteration=", iter, "\n")
    }
    
    # Check termination condition
    if ((dev_old - dev) / (0.1 + dev) < control[["dev_tol"]]) {
      if (control[["trace"]]) {
        cat("Convergence\n")
      }
      conv <- TRUE
      break
    }
  }
  
  # Compute Hessian, standard errors of \alpha, and null deviance
  H <- crossprod(MX)
  R <- try(chol(H), silent = TRUE)
  if (inherits(R, "try-error")) {
    se_alpha <- rep(Inf, ncol(H))
  } else {
    se_alpha <- sqrt(as.vector(variance_alpha(chol2inv(R), X, w_tilde^2, Ti)))
  }
  null_dev <- sum(family[["dev.resids"]](y, mean(y), wt))
  
  # Return list
  list(coefficients  = beta,
       alpha         = alpha,
       Hessian       = H,
       se_alpha      = se_alpha,
       deviance      = dev,
       null_deviance = null_dev,
       conv          = conv,
       iter          = iter)
}


# Offset function
bife_offset <- function(y, xb, id, Ti, family, control) {
  # Compute initial quantities for the maximization routine
  alpha <- family[["linkfun"]]((as.vector(tapply(y, id, mean)) + 0.5) / 2.0)
  eta <- alpha[id]
  alpha <- as.vector(tapply(eta - xb, id, mean))
  mu <- family[["linkinv"]](eta)
  wt <- rep(1.0, length(y))
  dev <- sum(family[["dev.resids"]](y, mu, wt))
  
  # Start maximization of the log-likelihood
  for (iter in seq.int(control[["iter_max"]])) {
    # Compute weights and dependent variable
    mu_eta <- family[["mu.eta"]](eta)
    w_tilde <- sqrt(mu_eta^2 / family[["variance"]](mu))
    nu <- (y - mu) / mu_eta
    
    # Compute update \alpha
    alpha_upd <- as.vector(update_alpha(nu * w_tilde, w_tilde, Ti))
    
    # Step-halving based on residual deviance as common in glm's
    dev_old <- dev
    rho <- 1.0
    repeat {
      # Compute residual deviance
      eta <- as.vector(xb + (alpha + rho * alpha_upd)[id])
      mu <- family[["linkinv"]](eta)
      dev <- sum(family[["dev.resids"]](y, mu, wt))
      
      # Check if deviance is not increasing
      if (is.finite(dev) && dev <= dev_old + control[["conv_tol"]] * dev) {
        # Update \alpha and leave step-halving
        alpha <- alpha + rho * alpha_upd
        break
      }
      
      # Update \rho
      rho <- rho / 2.0
      
      # If \rho is to small throw error
      if (rho < control[["rho_tol"]]) {
        stop("Failure in step-halving.", call. = FALSE)
      }
    }
    
    # Check termination condition
    if ((dev_old - dev) / (0.1 + dev) < control[["dev_tol"]]) {
      break
    }
  }
  
  # Return \alpha
  alpha
}
#' @title
#' Extract estimates of structural parameters or fixed effects
#' @description
#' \code{\link{coef.bife}} is a generic function which extracts estimates of the structural 
#' parameters or fixed effects from objects returned by \code{\link{bife}}.
#' @param 
#' object an object of class \code{"bife"}.
#' @param
#' type the type of parameter estimates that should be returned; structural parameters or 
#' fixed effects. Default is \code{"sp"} referring to the structural parameters.
#' @param 
#' corrected,fixed deprecated.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{\link{coef.bife}} returns a named vector of estimates of the requested 
#' parameters.
#' @seealso
#' \code{\link{bife}}
#' @export
coef.bife <- function(object, type = c("sp", "fe"), corrected = NULL, fixed = NULL, ...) {
  # Check validity of 'type'
  type <- match.arg(type)
  
  # 'corrected' is deprecated
  if (!is.null(corrected)) {
    warning(paste0("'corrected' is deprecated;",
                   "bias-corrected estimates are now automatically returned ", 
                   "if 'coef' is applied to an object returned by 'bias_corr'."), call. = FALSE)
  }
  
  # 'fixed' is deprecated
  if (!is.null(fixed)) {
    warning("'fixed' is deprecated; please use 'type = 'fe'' instead.", call. = FALSE)
    type <- "fe"
  }
  
  # Return requested estimates
  if (type == "sp") {
    object[["coefficients"]]
  } else {
    object[["alpha"]]
  }
}


#' @title
#' Extract estimates of average partial effects
#' @description
#' \code{\link{coef.bifeAPEs}} is a generic function which extracts estimates of the average partial 
#' effects from objects returned by \code{\link{get_APEs}}.
#' @param 
#' object an object of class \code{"APEs"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{\link{coef.bifeAPEs}} returns a named vector of estimates of the average 
#' partial effects.
#' @seealso
#' \code{\link{get_APEs}}
#' @export
coef.bifeAPEs <- function(object, ...) {
  object[["delta"]]
}


#' @title
#' Extract \code{bife} fitted values 
#' @description
#' \code{\link{fitted.bife}} is a generic function which extracts fitted values from an object 
#' returned by \code{\link{bife}}.
#' @param 
#' object an object of class \code{"bife"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{\link{fitted.bife}} returns a vector of fitted values.
#' @seealso
#' \code{\link{bife}}
#' @export
fitted.bife <- function(object, ...) {
  X <- model.matrix(object[["formula"]], object[["data"]], rhs = 1L)[, - 1L, drop = FALSE]
  x <- as.vector(X %*% object[["coefficients"]])
  x <- x + object[["alpha"]][as.integer(object[["data"]][[ncol(object[["data"]])]])]
  object[["family"]][["linkinv"]](x)
}


#' @title
#' Predict method for \code{bife} fits
#' @description
#' \code{\link{predict.bife}} is a generic function which obtains predictions from an object 
#' returned by \code{\link{bife}}.
#' @param 
#' object an object of class \code{"bife"}.
#' @param
#' type the type of prediction required. \code{"link"} is on the scale of the linear predictor
#' whereas \code{"response"} is on the scale of the response variable. Default is \code{"link"}.
#' @param 
#' X_new a regressor matrix for predictions. If not supplied predictions are based on the regressor
#' matrix returned by the object \code{\link{bife}}. See \code{Details}.
#' @param 
#' alpha_new a scalar or vector of fixed effects. If not supplied predictions are based on the
#' vector of fixed effects returned by \code{\link{bife}}. See \code{Details}.
#' @param 
#' corrected deprecated.
#' @param 
#' ... other arguments
#' @details 
#' The model frame returned by the object \code{\link{bife}} only includes individuals that 
#' were not dropped before the fitting process (due to perfect classification). The linear 
#' predictors of perfectly classified observations are equal to \code{- Inf} or \code{Inf} whereas 
#' the predicted probabilities are equal to their response. In-sample predictions are only based on 
#' non-perfectly classified observations.
#' 
#' If \code{alpha_new} is supplied as a scalar the linear predictor is computed using the same 
#' value of the fixed effect for each observation. If \code{alpha_new} is supplied as a vector it 
#' has to be of same length as the rows of the corresponding regressor matrix.
#' @return
#' The function \code{predict.bife} returns a vector of predictions.
#' @seealso
#' \code{\link{bife}}
#' @export
predict.bife <- function(object, type = c("link", "response"),
                         X_new = NULL, alpha_new = NULL, corrected = NULL, ...) {
  # 'corrected' is deprecated
  if (!is.null(corrected)) {
    warning(paste0("'corrected' is deprecated;",
                   "predictions are now automatically based on bias-corrected estimates ", 
                   "if 'predict' is applied to an object returned by 'bias_corr'."), call. = FALSE)
  }
  
  # Check validity of 'type'
  type <- match.arg(type)
  
  # 'X' provided for out-of-sample prediction?
  beta <- object[["coefficients"]]
  if (is.null(X_new)) {
    X <- model.matrix(object[["formula"]], object[["data"]], rhs = 1L)[, - 1L, drop = FALSE]
  } else {
    X <- as.matrix(X_new)
    if (ncol(X) != length(beta)) {
      stop("'X_new' of wrong dimension.", call. = FALSE)
    }
  }
  x <- as.vector(X %*% beta)
  
  # 'alpha' provided for out-of-sample prediction?
  if (is.null(alpha_new)) {
    if (is.null(X_new)) {
      x <- x + object[["alpha"]][as.integer(object[["data"]][[ncol(object[["data"]])]])]
    } else {
      x <- x + mean(object[["alpha"]])
    }
  } else {
    alpha <- as.vector(alpha_new)
    if (length(alpha) != 1L) {
      if (length(alpha) != length(x)) {
        stop("'alpha_new' of wrong dimension.", call. = FALSE)
      }
    }
    x <- x + alpha
  }
  
  # Compute requested type of prediction
  if (type == "response") {
    x <- object[["family"]][["linkinv"]](x)
  }
  
  # Return prediction
  x
}


#' @title
#' Print \code{bife}
#' @description
#' \code{\link{print.bife}} is a generic function which displays some minimal information from 
#' objects returned by \code{\link{bife}}.
#' @param 
#' x an object of class \code{"bife"}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{bife}}
#' @export
print.bife <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(x[["family"]][["family"]], " - ",
      x[["family"]][["link"]], " link\n\n", sep = "")
  print(x[["coefficients"]], digits = digits)
}


#' @title
#' Print \code{bifeAPEs}
#' @description
#' \code{\link{print.bifeAPEs}} is a generic function which displays some minimal information from 
#' objects returned by \code{\link{get_APEs}}.
#' @param 
#' x an object of class \code{"bifeAPEs"}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{get_APEs}}
#' @export
print.bifeAPEs <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print(x[["delta"]], digits = digits)
}


#' @title
#' Print \code{summary.bife}
#' @description
#' \code{\link{print.summary.bife}} is a generic function which displays summary statistics from 
#' objects returned by \code{\link{summary.bife}}.
#' @param 
#' x an object of class \code{"summary.bife"}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{bife}}
#' @export
print.summary.bife <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(x[["family"]][["family"]], " - ",
      x[["family"]][["link"]], " link\n\n", sep = "")
  print(x[["formula"]])
  cat("\nEstimates:\n")
  printCoefmat(x[["cm"]], P.values = TRUE, has.Pvalue = TRUE, digits = digits)
  cat("\nresidual deviance= ", round(x[["deviance"]], 2L), ",\n", sep = "")
  cat("null deviance= ", round(x[["null_deviance"]], 2L), ",\n", sep = "")
  cat("nT= ", x[["nobs"]][["nobs"]], ", N= ", x[["levels"]], "\n", sep = "")
  if (x[["nobs"]][["nobs_na"]] > 0L | x[["nobs"]][["nobs_pc"]] > 0L) {
    cat("\n")
    if (x[["nobs"]][["nobs_na"]] > 0L) {
      cat("(", x[["nobs"]][["nobs_na"]], "observation(s) deleted due to missingness )\n")
    }
    if (x[["nobs"]][["nobs_pc"]] > 0L) {
      cat("(", x[["nobs"]][["nobs_pc"]], "observation(s) deleted due to perfect classification )\n")
    }
  }
  cat("\nNumber of Fisher Scoring Iterations:", x[["iter"]], "\n")
  if (!is.null(x[["mean_alpha"]])) {
    cat("\nAverage individual fixed effect= ", round(x[["mean_alpha"]], 3L), "\n", sep = "")
  }
}


#' @title
#' Summarizing models of class \code{bife}
#' @description
#' Summary statistics for objects of class \code{"bife"}.
#' @param 
#' object an object of class \code{"bife"}.
#' @param
#' type the type of parameter estimates the summary statistics are related to: structural 
#' parameters or fixed effects. Default is \code{"sp"} referring to the structural parameters.
#' @param
#' corrected,fixed deprecated.
#' @param 
#' ... other arguments.
#' @return
#' Returns an object of class \code{"summary.bife"} which is a list of summary statistics of
#' \code{object}.
#' @seealso
#' \code{\link{bife}}
#' @export
summary.bife <- function(object, type = c("sp", "fe"), corrected = NULL, fixed = NULL, ...) {
  # Check validity of 'type'
  type <- match.arg(type)
  
  # 'corrected' is deprecated
  if (!is.null(corrected)) {
    warning(paste0("'corrected' is deprecated;",
                   "bias-corrected estimates are now automatically returned ", 
                   "if 'coef' is applied to an object returned by 'bias_corr'."), call. = FALSE)
  }
  
  # 'fixed' is deprecated
  if (!is.null(fixed)) {
    warning("'fixed' is deprecated; please use 'type = 'fe'' instead.", call. = FALSE)
    type <- "fe"
  }
  
  # Use estimates and standard errors requested
  if (type == "sp") {
    est <- object[["coefficients"]]
    se <- sqrt(diag(vcov(object)))
  } else {
    est <- object[["alpha"]]
    se <- object[["se_alpha"]]
  }
  
  # Compute coefficent matrix
  z <- est / se
  p <- 2.0 * pnorm(- abs(z))
  cm <- cbind(est, se, z, p)
  rownames(cm) <- names(est)
  colnames(cm) <- c("Estimate", "Std. error", "z value", "Pr(> |z|)")
  
  # Generate result list
  res <- list(cm            = cm, 
              deviance      = object[["deviance"]],
              null_deviance = object[["null_deviance"]],
              iter          = object[["iter"]],
              nobs          = object[["nobs"]],
              formula       = object[["formula"]],
              family        = object[["family"]],
              levels        = length(object[["alpha"]]))
  
  # Add average \alpha if estimates of the structural parameters are requested
  if (type == "sp") {
    res[["mean_alpha"]] <- mean(object[["alpha"]])
  }
  
  # Return list
  structure(res, class = "summary.bife")
}


#' @title
#' Print \code{summary.bifeAPEs}
#' @description
#' \code{\link{print.summary.bifeAPEs}} is a generic function which displays summary statistics from 
#' objects returned by \code{\link{summary.bifeAPEs}}.
#' @param 
#' x an object of class \code{"summary.bifeAPEs"}.
#' @param 
#' digits unsigned integer indicating the number of decimal places. Default is 
#' \code{max(3L, getOption("digits") - 3L)}.
#' @param 
#' ... other arguments.
#' @seealso
#' \code{\link{get_APEs}}
#' @export
print.summary.bifeAPEs <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Estimates:\n")
  printCoefmat(x, P.values = TRUE, has.Pvalue = TRUE, digits = digits)
}


#' @title
#' Summarizing models of class \code{bifeAPEs}
#' @description
#' Summary statistics for objects of class \code{"bifeAPEs"}.
#' @param 
#' object an object of class \code{"bifeAPEs"}.
#' @param 
#' ... other arguments.
#' @return
#' Returns an object of class \code{"summary.bifeAPEs"} which is a list of summary statistics of 
#' \code{object}.
#' @seealso
#' \code{\link{get_APEs}}
#' @export
summary.bifeAPEs <- function(object, ...) {
  # Compute coefficent matrix
  est <- object[["delta"]]
  se <- sqrt(diag(object[["vcov"]]))
  z <- est / se
  p <- 2.0 * pnorm(- abs(z))
  cm <- cbind(est, se, z, p)  
  rownames(cm) <- names(est)
  colnames(cm) <- c("Estimate", "Std. error", "z value", "Pr(> |z|)")
  
  # Return coefficient matrix
  structure(cm, class = "summary.bifeAPEs")
}


#' @title
#' Extract estimates of the covariance matrix
#' @description
#' \code{\link{vcov.bife}} computes an estimate of the covariance matrix of the estimator of the
#' structural parameters from objects returned by \code{\link{bife}}. The estimate is obtained
#' using the inverse of the negative Hessian after convergence.
#' @param 
#' object an object of class \code{"bife"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{\link{vcov.bife}} returns a named matrix of covariance estimates.
#' @seealso
#' \code{\link{bife}}
#' @export
vcov.bife <- function(object, ...) {
  # Check if the Hessian is invertible and compute its inverse
  R <- try(chol(object[["Hessian"]]), silent = TRUE)
  if (inherits(R, "try-error")) {
    matrix(Inf, length(coef(object)), length(coef(object)))
  } else {
    chol2inv(R)
  }
}


#' @title
#' Extract estimates of the covariance matrix
#' @description
#' \code{\link{vcov.bifeAPEs}} computes an estimate of the covariance matrix of the estimator of the
#' average partial parameters from objects returned by \code{\link{get_APEs}}.
#' @param 
#' object an object of class \code{"bifeAPEs"}.
#' @param 
#' ... other arguments.
#' @return
#' The function \code{\link{vcov.bifeAPEs}} returns a named matrix of covariance estimates.
#' @seealso
#' \code{\link{get_APEs}}
#' @export
vcov.bifeAPEs <- function(object, ...) {
  object[["vcov"]]
}
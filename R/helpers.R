### Helper functions (not exported)


# Checks if variable is a factor and transforms if necessary
check_factor <- function(x) {
  if (is.factor(x)) {
    droplevels(x)
  } else {
    factor(x)
  }
}


# Higher order partial derivatives of \mu with respect to \eta
partial_mu_eta <- function(eta, family, order) {
  f <- family[["mu.eta"]](eta)
  if (order == 2L) {
    if (family[["link"]] == "logit") {
      f * (1.0 - 2.0 * family[["linkinv"]](eta))
    } else {
      - eta * f
    }
  } else {
    if (family[["link"]] == "logit") {
      f * ((1.0 - 2.0 * family[["linkinv"]](eta))^2 - 2.0 * f)
    } else {
      (eta^2 - 1.0) * f
    }
  }
}


# Returns suitable name for a temporary variable
temp_var <- function(data) {
  repeat {
    tmpvar <- paste0(sample(letters, 5L, replace = TRUE), collapse = "")
    if (!(tmpvar %in% colnames(data))) {
      break
    }
  }
  tmpvar
}


# Unload
.onUnload <- function(libpath) {
  library.dynam.unload("bife", libpath)
}
\name{NEWS}
\title{bife news}
\section{Changes in version 0.7.2}{
  \itemize{
    \item Replaced NEWS.txt with NEWS.Rd.
    \item Improved the finite sample performance of bias corrections for the average partial effects in case of perfectly classified observations.
    \item The default option 'n_pop' in \code{get_APEs()} has been changed. Now the estimated covariance consists of the delta method part only, i.e. correction factor = 0.
  }
}
\section{Changes in version 0.7.1}{
  \itemize{
    \item Step-halving is now similar to \code{glm.fit2()}.
    \item Fixed an error in the covariance (influence function) of \code{get_APEs()}.
    \item Added reference for the bias correction in 'citation()'.
    \item Updated some references in the documentation and vignette.
    \item Fixed some typos in the documentation and vignette.
  }
}
\section{Changes in version 0.7}{
  \itemize{
    \item Added a generic function \code{logLik()} that returns the sum of the log-likelihood.
    \item 'X_new' in \code{predict()} now also permits data.frames.
    \item \code{bife()} now permits to expand functions with \code{poly()} and \code{bs()}.
    \item Added option 'sampling_fe' to \code{get_APEs()} to impose simplifying assumptions when estimating the covariance matrix.
    \item Fixed an error in \code{get_APEs()} related to the estimation of the covariance.
    \item Fixed a bug in the internal function used to estimate spectral densities.
  }
}
\section{Changes in version 0.6}{
  \itemize{
    \item Major source code cleaning and revision.
    \item Vignette renamed and revised.
    \item The asymptotic bias-correction is ported to a post-estimation routine \code{bias_corr()}.
    \item Replaced the bias-correction of Hahn and Newey (2004) with the one of Fernandez-Val (2009), i.e. refined bias expressions plus dynamic feedback possible.
    \item The post-estimation routine \code{apeff_bife()} was renamed to \code{get_APEs()} and also returns an estimate of the covariance using the delta method.
    \item Fixed a bug where the number of perfectly classified observations was not reported correctly. This bug also caused the estimates of the average partial effects to be unreliable for small time horizons.
    \item Some performance tweaks.
    \item Removed 'acs' data to reduce overall package size.
  }
}
\section{Changes in version 0.5}{
  \itemize{
    \item Changed the handling of \code{factor()} in the formula interface. This should now be more in line with what users expect.
    \item \code{fixed()} is removed since it is redundant now.
    \item Realigned \code{print.summary.bife()}.
  }
}
\section{Changes in version 0.4}{
  \itemize{
    \item Fixed a memory-access error in 'bife.cpp'.
  }
}
\section{Changes in version 0.3}{
  \itemize{
    \item Added a generic function \code{predict()}.
    \item Additional downgrades to ensure compability with g++ 4.4.x.
    \item Downgraded compiled code from C++11 to C++98. This ensures compability with older compilers.
    \item Removed jackknife bias-correction from \code{bife()} and \code{apeff_bife()}.
    \item \code{apeff_bife()} now uses the full sample instead of a sub-sample of indiviuals with a varying response.
    \item Changed default argument 'discrete' of \code{apeff_bife()} to 'NULL'.
    \item Changed default argument 'bias_corr' of \code{apeff_bife()} to '"ana"'.
    \item Fixed an issue where \code{vcov()} was not able to distinguish between corrected and uncorrected coefficients.
  }
}
\section{Changes in version 0.2}{
  \itemize{
    \item Added a function \code{fixed()} to model additional fixed-effects. See documentation for further details.
    \item Re-structured compiled code to optimize package size.
    \item Changed some function and argument names to common naming convention. For instance \code{apeff.bife(..., bias.corr = "ana")} is now \code{apeff_bife(..., bias_corr = "ana")}.
    \item Downgraded compiled code from C++11 to C++98. This ensures compability with older compilers.
    \item Added comments and renamed variables such that they are in line with the notation of the corresponding working paper.
    \item Fixed some typos in documentation.
    \item Fixed an issue where \code{bife()} was not able to fit a model with just one explanatory variable.
    \item Corrected computation of degrees of freedom.
  }
}
\section{Changes in version 0.1}{
  \itemize{
    \item Initial release on CRAN.
  }
}
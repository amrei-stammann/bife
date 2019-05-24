#' @title 
#' Female labor force participation
#' @description 
#' The sample was obtained from the "Panel Study of Income Dynamics" and
#' contains information about \eqn{N = 1461} women that were observed over
#' \eqn{T = 9} years.
#' @format A data frame with 13,149 rows:
#' \describe{
#'   \item{ID}{individual identifier}
#'   \item{LFP}{labor force participation}
#'   \item{KID1}{# of kids aged between 0 and 2}
#'   \item{KID2}{# of kids aged between 3 and 5}
#'   \item{KID3}{# of kids aged between 6 and 17}
#'   \item{INCH}{income husband}
#'   \item{AGE}{age of woman}
#'   \item{TIME}{time identifier}}
#' @references
#' Hyslop, D. (1999). "State Dependence, Serial Correlation and Heterogeneity in Intertemporal 
#' Labor Force Participation of Married Women". Econometrica 67(6), 1255-1294.
#' @references
#' Fernandez-Val, I. (2009). "Fixed effects estimation of structural parameters and marginal 
#' effects in panel probit models". Journal of Econometrics 150(1), 71-85.
#' @seealso
#' \code{\link{bife}}
"psid"
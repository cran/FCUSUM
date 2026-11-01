#' Long-run variance estimator (Bartlett kernel)
#'
#' @description
#' Computes the long-run variance of a residual series with the Bartlett
#' kernel and one of three bandwidth rules.
#'
#' @param u Numeric vector of residuals.
#' @param bw Bandwidth rule:
#'   \describe{
#'     \item{\code{"nw"}}{Newey-West (1994) deterministic rule
#'       \eqn{\lfloor 4 (T/100)^{2/9} \rfloor} (default).}
#'     \item{\code{"kurozumi"}}{Andrews (1991) AR(1) plug-in bounded as in
#'       Kurozumi (2002) with \eqn{k = 0.7}; better size with highly
#'       persistent errors, lower power.}
#'     \item{\code{"andrews"}}{Unbounded Andrews (1991) AR(1) plug-in
#'       (\eqn{\hat\rho} capped at 0.97).}
#'   }
#'
#' @return A list with the long-run variance (\code{omega2}) and the number
#'   of autocovariance lags used (\code{lags}).
#'
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent
#' covariance matrix estimation. \emph{Econometrica}, 59(3), 817-858.
#'
#' Kurozumi, E. (2002). Testing for stationarity with a break.
#' \emph{Journal of Econometrics}, 108(1), 63-99.
#'
#' Newey, W. K., & West, K. D. (1994). Automatic lag selection in covariance
#' matrix estimation. \emph{Review of Economic Studies}, 61(4), 631-653.
#'
#' @keywords internal
#' @noRd
lrv_bartlett <- function(u, bw = c("nw", "kurozumi", "andrews")) {
  bw <- match.arg(bw)
  u <- u - mean(u)
  n <- length(u)
  andrews_S <- function(r) {
    1.1447 * (4 * r^2 * n / ((1 + r)^2 * (1 - r)^2))^(1 / 3)
  }

  S <- switch(bw,
    nw = 4 * (n / 100)^(2 / 9),
    andrews = {
      rho <- sum(u[-1] * u[-n]) / sum(u[-n]^2)
      andrews_S(max(min(rho, 0.97), -0.97))
    },
    kurozumi = {
      rho <- sum(u[-1] * u[-n]) / sum(u[-n]^2)
      rho <- max(min(rho, 0.999), -0.999)
      min(andrews_S(rho), andrews_S(0.7))
    }
  )
  L <- min(floor(S), n - 2)

  s <- sum(u^2) / n
  if (L >= 1) {
    for (l in seq_len(L)) {
      s <- s + 2 * (1 - l / (L + 1)) * sum(u[(l + 1):n] * u[1:(n - l)]) / n
    }
  }
  list(omega2 = max(s, .Machine$double.eps), lags = L)
}

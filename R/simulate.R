#' Simulate Null Distributions for the Fourier CUSUM Test
#'
#' @description
#' Simulates, under the null of cointegration without structural change, the
#' distribution of one of the statistics used by \code{\link{fcum}}, with
#' exactly the same estimation steps (DOLS, Fourier terms, Bartlett long-run
#' variance). It reproduces the built-in tables and gives finite-sample
#' critical values for any configuration.
#'
#' @param n Integer. Sample size (before DOLS trimming).
#' @param p Integer. Number of I(1) regressors.
#' @param type Statistic to simulate:
#'   \describe{
#'     \item{\code{"plain"}}{Xiao-Phillips CUSUM without Fourier terms.}
#'     \item{\code{"supF"}}{sup-F pre-test statistic for the Fourier terms,
#'       maximised over the frequency grid defined by \code{kstar} and
#'       \code{freq}.}
#'     \item{\code{"fixed"}}{Fourier CUSUM at the fixed frequency \code{k}
#'       (no search).}
#'   }
#' @param k Numeric. Fourier frequency for \code{type = "fixed"}.
#' @param kstar,freq,step Frequency grid for \code{type = "supF"}; see
#'   \code{\link{fcum}}.
#' @param q Integer. Number of DOLS leads and lags.
#' @param bw Bandwidth rule for the long-run variance (see \code{\link{fcum}}).
#' @param nsim Integer. Number of Monte Carlo replications.
#' @param probs Probabilities for the returned quantiles.
#' @param seed Optional integer seed.
#'
#' @return A list with the simulated \code{quantiles} and all simulated
#'   statistics (\code{draws}).
#'
#' @details
#' Under the null, \eqn{y_t = u_t} with \eqn{u_t} iid N(0, 1) and \eqn{x_t} a
#' \eqn{p}-dimensional Gaussian random walk independent of \eqn{u_t}. The
#' statistics are invariant to the cointegrating vector, the intercept and
#' the scale of \eqn{u_t}.
#'
#' @examples
#' fcum_cv_sim(n = 100, p = 1, type = "plain", nsim = 200, seed = 1)$quantiles
#' fcum_cv_sim(n = 100, p = 1, type = "supF", kstar = 2, step = 0.05,
#'             nsim = 100, seed = 1)$quantiles
#'
#' @export
fcum_cv_sim <- function(n, p = 1, type = c("plain", "supF", "fixed"),
                        k = 1, kstar = 3, freq = c("fractional", "integer"),
                        step = 0.01, q = 1,
                        bw = c("nw", "kurozumi", "andrews"), nsim = 2000,
                        probs = c(0.50, 0.75, 0.80, 0.85, 0.90, 0.925,
                                  0.95, 0.975, 0.99),
                        seed = NULL) {
  type <- match.arg(type)
  freq <- match.arg(freq)
  bw <- match.arg(bw)
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  grid <- .freq_grid(kstar, freq, step)

  draws <- vapply(seq_len(nsim), function(i) {
    x <- .rw(n, p)
    y <- stats::rnorm(n)
    if (type == "fixed") {
      d <- .dols_design(y, x, q)
      return(.fcum_fit(d$y, d$Z, d$tt, n, k, bw)$statistic)
    }
    if (type == "plain") {
      d <- .dols_design(y, x, q)
      u <- stats::lm.fit(d$Z, d$y)$residuals
      return(max(abs(cumsum(u))) /
               (sqrt(lrv_bartlett(u, bw)$omega2) * sqrt(length(u))))
    }
    .fcum_core(y, x, q, kstar, freq, step, bw, grid = grid)$F
  }, numeric(1))

  list(quantiles = stats::quantile(draws, probs, names = TRUE), draws = draws)
}

#' Internal building blocks for the Fourier CUSUM test
#'
#' @description
#' Helpers shared by \code{fcum()} and the critical-value simulator:
#' DOLS design construction, lead/lag selection, the vectorised frequency
#' search and the CUSUM statistic itself.
#'
#' @keywords internal
#' @name fcum-internal
#' @noRd
NULL

# Prepare a numeric matrix of regressors and check inputs ----------------------
.as_regressors <- function(x, n) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop("All columns of 'x' must be numeric", call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (!is.numeric(x)) stop("'x' must be numeric", call. = FALSE)
  x <- as.matrix(x)                      # keeps n x p shape (old bug: flattened)
  storage.mode(x) <- "double"
  if (nrow(x) != n) {
    stop("'y' and 'x' must have the same number of observations", call. = FALSE)
  }
  if (anyNA(x)) stop("'x' contains missing values", call. = FALSE)
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  x
}

# DOLS design: constant, levels of x, q leads and lags of diff(x) ---------------
# 'trim' fixes the estimation sample (trim + 1):(n - trim) so that models with
# different q are estimated on the same observations when selecting q.
.dols_design <- function(y, x, q, trim = q) {
  n <- length(y)
  p <- ncol(x)
  dx <- rbind(matrix(NA_real_, 1, p), diff(x))
  idx <- (trim + 2):(n - trim)
  Z <- cbind(const = 1, x[idx, , drop = FALSE])
  if (q >= 0) {
    for (l in -q:q) {
      blk <- dx[idx + l, , drop = FALSE]
      colnames(blk) <- paste0("d", colnames(x), ifelse(l < 0, "_lag", "_lead"), abs(l))
      if (l == 0) colnames(blk) <- paste0("d", colnames(x))
      Z <- cbind(Z, blk)
    }
  }
  list(y = y[idx], Z = Z, tt = idx)
}

# Select the number of DOLS leads/lags by BIC (common sample) -------------------
.select_q <- function(y, x, qmax) {
  bic <- vapply(0:qmax, function(q) {
    d <- .dols_design(y, x, q, trim = qmax)
    e <- stats::lm.fit(d$Z, d$y)$residuals
    m <- length(e)
    m * log(sum(e^2) / m) + ncol(d$Z) * log(m)
  }, numeric(1))
  (0:qmax)[which.min(bic)]
}

# Vectorised frequency search via Frisch-Waugh-Lovell ----------------------------
# For every candidate k, RSS of y ~ Z + cos(2 pi k t / n) + sin(2 pi k t / n).
# All candidates have the same number of parameters, so minimum RSS equals
# minimum AIC / AICc / BIC.
.freq_search <- function(yy, Z, tt, n, grid) {
  qrZ <- qr(Z)
  ey <- qr.resid(qrZ, yy)
  arg <- 2 * pi * outer(tt, grid) / n
  RC <- qr.resid(qrZ, cos(arg))
  RS <- qr.resid(qrZ, sin(arg))
  a <- colSums(RC^2)
  b <- colSums(RC * RS)
  d <- colSums(RS^2)
  g1 <- colSums(RC * ey)
  g2 <- colSums(RS * ey)
  det <- a * d - b^2
  ess <- (d * g1^2 - 2 * b * g1 * g2 + a * g2^2) / det
  rss <- sum(ey^2) - ess
  rss[!is.finite(rss) | det <= 1e-10 * a * d] <- Inf
  rss
}

# CUSUM statistic for a given frequency -------------------------------------------
.fcum_fit <- function(yy, Z, tt, n, k, bw) {
  X <- cbind(Z,
             cos_k = cos(2 * pi * k * tt / n),
             sin_k = sin(2 * pi * k * tt / n))
  fit <- stats::lm.fit(X, yy)
  u <- fit$residuals
  m <- length(u)
  lv <- lrv_bartlett(u, bw = bw)
  stat <- max(abs(cumsum(u))) / (sqrt(lv$omega2) * sqrt(m))
  list(statistic = stat, residuals = u, coefficients = fit$coefficients,
       X = X, omega2 = lv$omega2, bw_lags = lv$lags, nobs = m)
}

# Full procedure on one sample -----------------------------------------------------
# DOLS without and with the Fourier pair at the RSS-minimising frequency.
# Returns both CUSUM statistics and the sup-F statistic for the Fourier terms:
#   F = [(RSS_0 - RSS_k)/2] / omega2_k,
# where omega2_k is the long-run variance of the Fourier-model residuals.
# Because all Fourier candidates share RSS_0 and have 2 extra parameters, the
# RSS-minimising frequency also maximises the numerator of F.
.fcum_core <- function(y, x, q, kstar, freq, step, bw, grid = NULL) {
  n <- length(y)
  d <- .dols_design(y, x, q)
  if (is.null(grid)) grid <- .freq_grid(kstar, freq, step)
  rss <- .freq_search(d$y, d$Z, d$tt, n, grid)
  khat <- grid[which.min(rss)]

  f1 <- .fcum_fit(d$y, d$Z, d$tt, n, khat, bw)            # with Fourier pair
  f0 <- stats::lm.fit(d$Z, d$y)                           # no Fourier terms
  u0 <- f0$residuals
  m <- length(u0)
  lv0 <- lrv_bartlett(u0, bw = bw)
  stat0 <- max(abs(cumsum(u0))) / (sqrt(lv0$omega2) * sqrt(m))
  Fstat <- ((sum(u0^2) - sum(f1$residuals^2)) / 2) / f1$omega2

  list(khat = khat, fourier = f1, stat_fourier = f1$statistic,
       stat_plain = stat0, F = Fstat, plain_coef = f0$coefficients,
       plain_resid = u0, plain_omega2 = lv0$omega2, plain_bw_lags = lv0$lags,
       n = n, nobs = m, grid = grid)
}

# Simulate p independent Gaussian random walks -----------------------------------------
.rw <- function(n, p) {
  x <- apply(matrix(stats::rnorm(n * p), n, p), 2, cumsum)
  x <- matrix(x, n, p)
  colnames(x) <- paste0("x", seq_len(p))
  x
}

# Bootstrap null distribution of the Fourier CUSUM statistic ----------------------------
# Data are generated under H0 (cointegration) with a Fourier break at the
# estimated frequency and with the estimated standardised Fourier amplitudes,
# and the whole procedure (DOLS, frequency search, LRV) is re-run on each
# bootstrap sample. The statistic is invariant to the cointegrating vector and
# to the scale of u, so x* is a standard random walk and u* is standardised.
#   type = "iid":   u* iid N(0, 1)
#   type = "sieve": u* from an AR(p) (BIC, p <= 4) fitted to the residuals,
#                   with resampled innovations; mimics short-run dynamics.
.boot_fourier <- function(core, p, q, bw, nboot, type) {
  n <- core$n
  tt <- seq_len(n)
  fit <- core$fourier
  sc <- sqrt(fit$omega2)
  amp <- unname(utils::tail(fit$coefficients, 2)) / sc
  k <- core$khat
  brk <- amp[1] * cos(2 * pi * k * tt / n) + amp[2] * sin(2 * pi * k * tt / n)

  if (type == "sieve") {
    ar <- stats::ar(fit$residuals, aic = TRUE, order.max = 4, method = "ols",
                    demean = TRUE)
    phi <- as.numeric(ar$ar)
    if (length(phi)) {
      s <- sum(phi)
      if (s > 0.9) phi <- phi * 0.9 / s                 # keep it stationary
      if (any(Mod(polyroot(c(1, -phi))) < 1.02)) phi <- numeric(0)
    }
    e <- as.numeric(stats::na.omit(ar$resid))
    e <- (e - mean(e)) / sc
  }

  vapply(seq_len(nboot), function(b) {
    x <- .rw(n, p)
    if (type == "iid") {
      u <- stats::rnorm(n)
    } else {
      burn <- 50
      e_b <- sample(e, n + burn, replace = TRUE)
      u <- if (length(phi)) {
        as.numeric(stats::filter(e_b, phi, method = "recursive"))
      } else e_b
      u <- utils::tail(u, n)
    }
    y <- brk + u
    d <- .dols_design(y, x, q)
    rss <- .freq_search(d$y, d$Z, d$tt, n, core$grid)
    kb <- core$grid[which.min(rss)]
    .fcum_fit(d$y, d$Z, d$tt, n, kb, bw)$statistic
  }, numeric(1))
}

# Frequency grid ------------------------------------------------------------------
.freq_grid <- function(kstar, freq, step) {
  if (freq == "integer") {
    seq_len(kstar)
  } else {
    round(seq(0.1, kstar, by = step), 10)
  }
}

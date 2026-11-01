#' Print Method for fcum Objects
#'
#' @param x An object of class \code{fcum}.
#' @param digits Number of digits to print.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.fcum <- function(x, digits = 4, ...) {
  f <- function(v) format(round(v, digits), nsmall = digits)
  s <- x$settings
  cat("\nFourier CUSUM Test for the Null of Cointegration\n")
  cat("=================================================\n")
  cat("H0: cointegration  |  H1: no cointegration\n\n")

  pt <- x$pretest
  cat("Step 1 - Fourier pre-test (sup-F)\n")
  cat("  sup-F = ", f(pt$statistic), "  p-value: ", pt$p.label,
      "  (", 100 * pt$level, "% c.v. = ", f(pt$cv_level), ")\n", sep = "")
  cat("  Selected frequency k-hat: ", x$best_frequency, " (", s$freq,
      " grid, k* = ", s$kstar, ")\n", sep = "")
  cat("  Fourier terms ",
      if (pt$fourier_retained) "retained" else "dropped",
      if (pt$mode != "auto") paste0(" (fourier = \"", pt$mode, "\")") else "",
      "\n\n", sep = "")

  cat("Step 2 - ", x$test, " test\n", sep = "")
  cat("  CUSUM statistic: ", f(x$statistic), " ", x$significance, "\n", sep = "")
  cat("  p-value:         ", x$p.label, "\n", sep = "")
  cat("  Critical values (", x$cv_source, "):\n", sep = "")
  cat("    10%: ", f(x$critical_values["10%"]),
      "   5%: ", f(x$critical_values["5%"]),
      "   1%: ", f(x$critical_values["1%"]), "\n\n", sep = "")
  cat("Decision:", x$decision, "\n\n")
  cat("Regressors: ", s$p, "  |  DOLS leads/lags: ", s$q,
      "  |  Obs. used: ", x$nobs, "  |  Bandwidth: ", s$bw, "\n", sep = "")
  invisible(x)
}

#' Summary Method for fcum Objects
#'
#' @param object An object of class \code{fcum}.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#'
#' @export
summary.fcum <- function(object, ...) {
  print(object)
  cat("\nCall:\n")
  print(object$call)
  cat("\nBoth CUSUM statistics:  Fourier = ", round(object$stat_fourier, 4),
      "   Xiao-Phillips (no Fourier) = ", round(object$stat_plain, 4), "\n",
      sep = "")
  cat("Long-run variance of residuals: ", signif(object$lrv, 5),
      "  (Bartlett, ", object$settings$bw, ", ", object$bw_lags, " lags)\n",
      sep = "")
  cat("\nDOLS coefficients (", object$test, " model):\n", sep = "")
  co <- object$coefficients
  print(round(co[!grepl("^d", names(co))], 4))
  cat("(lead/lag terms of diff(x) not shown; the default lm standard errors\n",
      " in object$best_model are not valid for inference with DOLS)\n", sep = "")
  invisible(object)
}

#' Plot Method for fcum Objects
#'
#' @description
#' Plots the scaled cumulative sum of residuals of the reported model with
#' the critical band.
#'
#' @param x An object of class \code{fcum}.
#' @param level Significance level for the band: \code{"5\%"} (default),
#'   \code{"10\%"} or \code{"1\%"}.
#' @param ... Further graphical parameters passed to \code{plot}.
#'
#' @return Invisibly returns the input object.
#'
#' @export
plot.fcum <- function(x, level = c("5%", "10%", "1%"), ...) {
  level <- match.arg(level)
  path <- cumsum(x$residuals) / (sqrt(x$lrv) * sqrt(x$nobs))
  cvl <- x$critical_values[level]
  lim <- range(c(path, cvl, -cvl))
  graphics::plot(seq_along(path), path, type = "l", ylim = lim,
                 xlab = "Observation", ylab = "Scaled CUSUM of residuals",
                 main = paste(x$test, "path"), ...)
  graphics::abline(h = c(-cvl, cvl), lty = 2, col = "red")
  graphics::abline(h = 0, col = "grey")
  graphics::legend("topleft", bty = "n", lty = 2, col = "red",
                   legend = paste0(level, " critical value"))
  invisible(x)
}

#' @rdname summary.fcum
#' @export
coef.fcum <- function(object, ...) object$coefficients

#' @rdname summary.fcum
#' @export
residuals.fcum <- function(object, ...) object$residuals

#' Fourier CUSUM Test for the Null of Cointegration
#'
#' @description
#' Tests the null hypothesis of cointegration against the alternative of no
#' cointegration, allowing for smooth structural change in the intercept
#' approximated by a Fourier pair at a fractional (default) or integer
#' frequency. The procedure extends the residual-based CUSUM test of Xiao
#' and Phillips (2002) with the flexible Fourier form, in three steps:
#' \enumerate{
#'   \item estimate the cointegrating regression by dynamic OLS and search the
#'     Fourier frequency that minimises the residual sum of squares;
#'   \item pre-test the Fourier terms with a sup-F statistic;
#'   \item if the Fourier terms are significant, compute the Fourier CUSUM
#'     statistic with bootstrap critical values that replicate the frequency
#'     search; otherwise compute the Xiao-Phillips CUSUM statistic without
#'     Fourier terms.
#' }
#'
#' @param y Numeric vector. Dependent I(1) variable.
#' @param x Numeric vector, matrix or data frame of I(1) regressors with the
#'   same number of observations as \code{y}.
#' @param kstar Positive integer. Maximum Fourier frequency searched
#'   (default 3).
#' @param freq \code{"fractional"} (default; grid \code{0.1, 0.1 + step,
#'   ..., kstar}) or \code{"integer"} (grid \code{1, ..., kstar}).
#' @param fourier \code{"auto"} (default) decides with the sup-F pre-test;
#'   \code{"always"} reports the Fourier CUSUM regardless; \code{"never"}
#'   reports the Xiao-Phillips CUSUM without Fourier terms.
#' @param pretest_level Significance level of the sup-F pre-test
#'   (default 0.05).
#' @param cv Critical values for the Fourier CUSUM: \code{"bootstrap"}
#'   (default) or \code{"table"} (asymptotic values at the selected
#'   frequency; fast but conservative in finite samples because they ignore
#'   the frequency search).
#' @param boot Bootstrap scheme: \code{"iid"} (default; best power) or
#'   \code{"sieve"} (AR sieve on the residuals; better size with highly
#'   persistent errors, lower power).
#' @param nboot Integer. Bootstrap replications (default 499). Also used
#'   for simulated critical values when a configuration is outside the
#'   built-in tables.
#' @param q Integer or \code{NULL}. DOLS leads and lags of \eqn{\Delta x_t};
#'   chosen by BIC between 0 and \code{qmax} when \code{NULL} (default).
#' @param qmax Integer. Largest number of leads/lags considered. Default
#'   \code{min(4, floor(4 * (n / 100)^(1/4)))}.
#' @param bw Bandwidth rule for the Bartlett long-run variance:
#'   \code{"nw"} (Newey-West, default), \code{"kurozumi"} (Andrews AR(1)
#'   plug-in bounded as in Kurozumi 2002) or \code{"andrews"}.
#' @param step Grid step for fractional frequencies (default 0.01).
#' @param seed Optional integer seed for the bootstrap.
#'
#' @return An object of class \code{fcum}, a list with
#' \describe{
#'   \item{test}{\code{"Fourier CUSUM"} or \code{"Xiao-Phillips CUSUM"}:
#'     the statistic on which the decision is based.}
#'   \item{statistic, p.value, critical_values, decision, significance}{For
#'     the reported test.}
#'   \item{pretest}{List with the sup-F statistic, its critical values,
#'     p-value and whether the Fourier terms were retained.}
#'   \item{stat_fourier, stat_plain}{Both CUSUM statistics.}
#'   \item{best_frequency}{Selected Fourier frequency.}
#'   \item{coefficients, residuals, lrv, bw_lags, nobs}{Of the model behind
#'     the reported test.}
#'   \item{best_model}{That model as an \code{lm} object (its default
#'     standard errors are not valid for DOLS inference).}
#'   \item{settings}{kstar, freq, q, bw, cv, boot, nboot, p.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details
#' The Fourier cointegrating regression is
#' \deqn{y_t = \alpha_0 + \alpha_1 \sin(2\pi k t/T) + \alpha_2 \cos(2\pi k t/T)
#'       + \beta' x_t + \sum_{j=-q}^{q} \gamma_j' \Delta x_{t-j} + u_t,}
#' estimated by dynamic OLS so that the limit distribution of the statistic
#' is free of nuisance parameters arising from endogeneity and serial
#' correlation (Xiao and Phillips 2002). The CUSUM statistic is
#' \deqn{CS = \max_{t} \frac{1}{\hat{\omega}\sqrt{T}}
#'       \left| \sum_{s=1}^{t} \hat{u}_s \right|,}
#' with \eqn{\hat\omega^2} a Bartlett long-run variance estimate. Large values
#' reject cointegration.
#'
#' \strong{Why a pre-test and a bootstrap.} With a fractional frequency the
#' null distribution of the Fourier CUSUM depends on the true frequency, and
#' the search over \eqn{k} shrinks the statistic. When there is no break the
#' selected frequency is spurious; when there is one, it is estimated with
#' error. The sup-F pre-test (Becker, Enders and Lee 2006; Tsong et al. 2016)
#' sends samples without a break to the standard Xiao-Phillips test, and the
#' parametric bootstrap reproduces the frequency search around the estimated
#' break so that the Fourier CUSUM has correct size. In Monte Carlo
#' experiments (see \code{data-raw/size_power.R}) this keeps size close to the
#' nominal level both with and without a break, whereas critical values
#' conditional on the selected frequency alone are very conservative and
#' have little power.
#'
#' @section Changes from version 1.0.0:
#' The 1.0.0 implementation (i) used critical values that were lower-tail
#' quantiles, so a true null was rejected in about 90\% of samples at the
#' nominal 5\% level; (ii) used OLS residuals scaled by their short-run
#' standard deviation, so the test was not robust to serial correlation or
#' endogeneity; (iii) used critical values indexed by \code{kstar} that
#' ignore the dependence on the frequency; and (iv) flattened a multivariate
#' \code{x} into one column. All are corrected here.
#'
#' @references
#' Xiao, Z., & Phillips, P. C. B. (2002). A CUSUM test for cointegration
#' using regression residuals. \emph{Journal of Econometrics}, 108(1), 43-61.
#'
#' Tsong, C.-C., Lee, C.-F., Tsai, L.-J., & Hu, T.-C. (2016). The Fourier
#' approximation and testing for the null of cointegration.
#' \emph{Empirical Economics}, 51(3), 1085-1113.
#'
#' Becker, R., Enders, W., & Lee, J. (2006). A stationarity test in the
#' presence of an unknown number of smooth breaks. \emph{Journal of Time
#' Series Analysis}, 27(3), 381-409.
#'
#' Stock, J. H., & Watson, M. W. (1993). A simple estimator of cointegrating
#' vectors in higher order integrated systems. \emph{Econometrica}, 61(4),
#' 783-820.
#'
#' Kurozumi, E. (2002). Testing for stationarity with a break.
#' \emph{Journal of Econometrics}, 108(1), 63-99.
#'
#' @examples
#' set.seed(123)
#' n <- 150
#' x <- cumsum(rnorm(n))
#' brk <- 2 * sin(2 * pi * 1.5 * seq_len(n) / n)
#' y <- 2 + brk + 1.5 * x + rnorm(n)
#'
#' res <- fcum(y, x, kstar = 3, nboot = 99, seed = 1)
#' res
#' summary(res)
#'
#' # Two regressors, integer frequencies, fast table critical values
#' x2 <- cbind(x1 = x, x2 = cumsum(rnorm(n)))
#' y2 <- 1 + brk + x2 %*% c(1, -0.5) + rnorm(n)
#' fcum(y2, x2, kstar = 2, freq = "integer", cv = "table")
#'
#' @export
fcum <- function(y, x, kstar = 3, freq = c("fractional", "integer"),
                 fourier = c("auto", "always", "never"), pretest_level = 0.05,
                 cv = c("bootstrap", "table"), boot = c("iid", "sieve"),
                 nboot = 499, q = NULL, qmax = NULL,
                 bw = c("nw", "kurozumi", "andrews"), step = 0.01,
                 seed = NULL) {
  freq <- match.arg(freq)
  fourier <- match.arg(fourier)
  cv <- match.arg(cv)
  boot <- match.arg(boot)
  bw <- match.arg(bw)

  # ---- input checks ----------------------------------------------------------
  if (missing(y) || missing(x)) stop("Arguments 'y' and 'x' are required", call. = FALSE)
  if (is.data.frame(y) || is.matrix(y)) {
    if (NCOL(y) != 1) stop("'y' must be a single series", call. = FALSE)
    y <- y[, 1]
  }
  if (!is.numeric(y)) stop("'y' must be numeric", call. = FALSE)
  y <- as.numeric(y)
  if (anyNA(y)) stop("'y' contains missing values", call. = FALSE)
  n <- length(y)
  x <- .as_regressors(x, n)
  p <- ncol(x)

  if (!is.numeric(kstar) || length(kstar) != 1 || kstar < 1 ||
      kstar != round(kstar)) {
    stop("'kstar' must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(step) || step <= 0 || step > 1) {
    stop("'step' must be in (0, 1]", call. = FALSE)
  }
  if (!is.numeric(pretest_level) || pretest_level <= 0 || pretest_level >= 1) {
    stop("'pretest_level' must be in (0, 1)", call. = FALSE)
  }
  if (!is.numeric(nboot) || nboot < 19) stop("'nboot' must be at least 19", call. = FALSE)
  if (is.null(qmax)) qmax <- min(4L, floor(4 * (n / 100)^(1 / 4)))
  if (!is.null(q) && (q < 0 || q != round(q))) {
    stop("'q' must be a non-negative integer", call. = FALSE)
  }

  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }

  # ---- step 1: DOLS + frequency search -------------------------------------------
  if (is.null(q)) q <- .select_q(y, x, qmax)
  d <- .dols_design(y, x, q)
  if (length(d$y) - ncol(d$Z) - 2 < 10) {
    stop("Too few observations for ", p, " regressor(s) and q = ", q,
         " leads/lags", call. = FALSE)
  }
  core <- .fcum_core(y, x, q, kstar, freq, step, bw)

  # ---- step 2: sup-F pre-test ------------------------------------------------------
  tabF <- if (step == 0.01 || freq == "integer") cv_supF(p, kstar, freq) else NULL
  if (is.null(tabF)) {
    message("sup-F critical values not tabulated for this configuration; ",
            "simulating (nsim = ", nboot, ").")
    simF <- fcum_cv_sim(n, p, type = "supF", kstar = kstar, freq = freq,
                        step = step, q = q, bw = bw, nsim = nboot)
    csF <- cv_summary(core$F, draws = simF$draws)
  } else {
    csF <- cv_summary(core$F, tabF$quantiles, tabF$probs)
  }
  cvF_level <- if (is.null(tabF)) {
    stats::quantile(simF$draws, 1 - pretest_level, names = FALSE)
  } else {
    stats::approx(tabF$probs, tabF$quantiles, xout = 1 - pretest_level,
                  rule = 2)$y
  }
  use_fourier <- switch(fourier,
                        auto = core$F > cvF_level,
                        always = TRUE,
                        never = FALSE)
  pretest <- list(statistic = c(supF = core$F),
                  critical_values = csF$critical_values,
                  p.value = csF$p.value, p.label = csF$p.label,
                  level = pretest_level, cv_level = cvF_level,
                  fourier_retained = use_fourier, mode = fourier)

  # ---- step 3: CUSUM test ---------------------------------------------------------
  if (use_fourier) {
    test <- "Fourier CUSUM"
    stat <- core$stat_fourier
    if (cv == "bootstrap") {
      draws <- .boot_fourier(core, p, q, bw, nboot, boot)
      cs <- cv_summary(stat, draws = draws)
      cv_source <- paste0(boot, " bootstrap (", nboot, " reps)")
    } else {
      tq <- cv_fixed(p, core$khat)
      if (is.null(tq)) {
        message("Outside the fixed-frequency table; simulating (nsim = ", nboot, ").")
        sim <- fcum_cv_sim(n, p, type = "fixed", k = core$khat, q = q, bw = bw,
                           nsim = nboot)
        cs <- cv_summary(stat, draws = sim$draws)
        cv_source <- "simulated, fixed frequency"
      } else {
        cs <- cv_summary(stat, tq$quantiles, tq$probs)
        cv_source <- "asymptotic table, fixed frequency"
      }
    }
    fit <- core$fourier
    coefs <- fit$coefficients; resid <- fit$residuals
    lrv <- fit$omega2; bw_lags <- fit$bw_lags; X <- fit$X
  } else {
    test <- "Xiao-Phillips CUSUM"
    stat <- core$stat_plain
    tq <- cv_plain(p)
    if (is.null(tq)) {
      message("Plain CUSUM table covers p <= 4; simulating (nsim = ", nboot, ").")
      sim <- fcum_cv_sim(n, p, type = "plain", q = q, bw = bw, nsim = nboot)
      cs <- cv_summary(stat, draws = sim$draws)
      cv_source <- "simulated"
    } else {
      cs <- cv_summary(stat, tq$quantiles, tq$probs)
      cv_source <- "asymptotic table"
    }
    coefs <- core$plain_coef; resid <- core$plain_resid
    lrv <- core$plain_omega2; bw_lags <- core$plain_bw_lags; X <- d$Z
  }
  dec <- .decide(stat, cs$critical_values)

  # ---- lm object for display -------------------------------------------------------
  Xd <- as.data.frame(X[, -1, drop = FALSE])
  names(Xd) <- make.names(names(Xd), unique = TRUE)
  Xd$.y <- d$y
  best_model <- stats::lm(.y ~ ., data = Xd)

  structure(
    list(
      test = test,
      statistic = c(CUSUM = unname(stat)),
      p.value = cs$p.value,
      p.label = cs$p.label,
      critical_values = cs$critical_values,
      cv_source = cv_source,
      decision = dec$decision,
      significance = dec$significance,
      pretest = pretest,
      stat_fourier = unname(core$stat_fourier),
      stat_plain = unname(core$stat_plain),
      best_frequency = core$khat,
      coefficients = coefs,
      residuals = resid,
      lrv = lrv,
      bw_lags = bw_lags,
      nobs = core$nobs,
      p = p,
      best_model = best_model,
      settings = list(kstar = kstar, freq = freq, q = q, bw = bw, cv = cv,
                      boot = boot, nboot = nboot, step = step, p = p),
      call = match.call()
    ),
    class = "fcum"
  )
}

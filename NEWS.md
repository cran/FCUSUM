# FCUSUM 1.1.0

## Critical corrections — results from 1.0.0 should be recomputed

* **Critical values were the wrong tail.** The table in 1.0.0 contained the
  1%, 5% and 10% *lower* quantiles of the null distribution. Because the
  test rejects for large values, `fcum()` rejected a true null of
  cointegration in roughly 90% of samples at the nominal 5% level.
* **Estimation.** The cointegrating regression is now estimated by dynamic
  OLS (leads and lags of the differenced regressors; number chosen by BIC
  unless `q` is given), as Xiao and Phillips (2002) require for a
  nuisance-parameter-free limit distribution. 1.0.0 used plain OLS.
* **Scaling.** The CUSUM is now scaled by a Bartlett long-run variance
  (`bw = "nw"` default, `"kurozumi"`, `"andrews"`) instead of the residual
  standard deviation, so the test is robust to serially correlated errors.
* **Multiple regressors.** `x` given as a matrix or data frame was
  flattened into a single column in 1.0.0, so p > 1 failed. Fixed.

## New testing procedure

With fractional frequencies the null distribution of the Fourier CUSUM
depends on the true frequency, and the frequency search shrinks the
statistic, so no single table indexed by `kstar` can be valid. `fcum()` now
proceeds in three steps:

1. DOLS estimation and search for the frequency that minimises the RSS.
2. A sup-F pre-test for the Fourier terms (tabulated critical values for
   p = 1..4, kstar = 1..5, fractional and integer grids).
3. If the Fourier terms are significant, the Fourier CUSUM with bootstrap
   critical values that replicate the frequency search around the estimated
   break (`boot = "iid"`, default, or `"sieve"`); otherwise the
   Xiao-Phillips CUSUM without Fourier terms (tabulated critical values).

`fourier = "always"` / `"never"` override the pre-test; `cv = "table"` uses
fast asymptotic critical values at the selected frequency instead of the
bootstrap (conservative in finite samples). Configurations outside the
tables are simulated automatically.

## Other changes

* `kstar` must be a positive integer (non-integer values silently fell back
  to the p = 1, k = 1 critical values in 1.0.0).
* Frequency selection is by minimum RSS. All candidate models have the same
  number of parameters, so this is identical to the AICc rule of 1.0.0; the
  internal `aicc()` helper was removed.
* New exported function `fcum_cv_sim()` simulates the null distribution of
  the plain CUSUM, the sup-F statistic or the fixed-frequency Fourier CUSUM
  for any n, p and settings.
* `fcum` objects now carry p-values, the pre-test, both CUSUM statistics,
  and gain `plot()`, `coef()` and `residuals()` methods.
* `data-raw/critical_values.R` reproduces all tables and
  `data-raw/size_power.R` reproduces the size and power experiments.

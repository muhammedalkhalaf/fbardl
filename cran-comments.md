## fbardl 1.2.0

This release corrects the computations below; the version on CRAN is 1.1.0.

It follows 1.1.0 within days because a re-audit found that the bootstrap held the Fourier frequency and the lag selection fixed at their data-selected values (which made the bootstrap tests oversized) and that the analytical bounds verdict was invalid with Fourier terms, so the results users obtain change.

## Breaking changes

* The default `type` is now `"fbardl_bvz"` (was `"fardl"`). With the default `fourier = TRUE` the Kripfganz and Schneider bounds do not account for the Fourier terms (Monte Carlo size of the 5% bounds decision 0.605 with independent random walks, n = 100), so the default call now returns a bootstrap test. The Bertelli, Vacca and Zoia scheme is the default for consistency with ardlverse::fbnardl. Use `type = "fardl"` for the previous behaviour.
* The default grid of candidate Fourier frequencies is now `kgrid = "integer"`, k = 1, ..., `maxk` (the integer frequencies selected by minimum SSR in Enders and Lee, 2012); before it was the fractional grid 0.1, 0.2, ..., `maxk`, which is still available as `kgrid = "fractional"`. With the default arguments the selected k* can change, and with it every estimate and statistic. `kgrid = "fractional"` reproduces fbardl 1.1.0 exactly. The integer default, and whether the fractional grid is appropriate, are package choices pending verification against Yilanci, Bozoklu and Gorus (2020) and Omay (2015).
* `type = "fardl"` with Fourier terms no longer gives a decision or p-values: the bounds are shown labelled "bounds for the model without Fourier terms; not valid with Fourier terms", and the number of short-run coefficients used for them now excludes the two Fourier terms (so the displayed critical values change). `print()` and `summary()` say so and point to the bootstrap types. Without Fourier terms (`fourier = FALSE`) the bounds, p-values and decision are unchanged.
* Bootstrap decisions use the labels of the shared engine: "Cointegration" requires Fov, t and Find to reject at 5%; "degenerate case of the first type" (Fov rejects and t does not, whatever Find gives) and "of the second type" (Fov and t reject, Find does not), in the terminology of Bertelli, Vacca and Zoia (2022, eqs. 8-9); McNown, Sam and Goh (2018) number the two degenerate cases the other way. The same labels are used for both bootstrap types; otherwise "No cointegration". The 1.1.0 label "PARTIAL EVIDENCE" is gone. A code is returned in `cointegration$decision.code`.

## Bootstrap corrections

* The bootstrap types now use the recursive bootstrap engine shared with ardlverse (file R/ardl_boot_engine.R, engine version 1.1.0, copied verbatim and checked by a test). fbardl 1.1.0 selected k* and the lags on the data and held them fixed in every replication, which made the bootstrap tests oversized for the Fourier model. Now, in every replication, k* is selected again on the bootstrap sample (unless fixed with the new argument `kstar`, or `fourier = FALSE`) and the lags are selected again (unless fixed with the new argument `lags`); the model is then estimated with the selected specification and Fov, t and Find are computed exactly as on the data. Repeating the selection is a package choice; BVZ step 5(c) re-estimates the unrestricted model but does not discuss re-selection, and McNown, Sam and Goh do not discuss it either.
* `"fbardl_bvz"`: separate nulls for Fov, t and Find, marginal VECM for the differences of the regressors, residuals recentred after each draw (Bertelli, Vacca and Zoia, 2022, Section 3, as before) and, new, initial values drawn as a random block of the data (their step 5b); 1.1.0 used the first observations.
* `"fbardl_mcnown"`: the Fov null for all statistics (McNown, Sam and Goh, 2018, Steps 1 to 8), applied to the conditional ECM by default (a package choice; MSG write the y equation in the unconditional form of their eq. 12, used with `unconditional = TRUE`), with the unrestricted equation for the differences of the regressors and observed initial values; residuals are recentred once by their mean (MSG eq. 13; 1.1.0 divided their sum by the number of residuals minus p minus 1, following the Stata module).
* Critical values are the order statistics of MSG eqs. 15-16 and BVZ eqs. 24-25 (the lower t critical value is now the (floor(0.05 B) + 1)-th smallest statistic; 1.1.0 took the floor(0.05 B)-th). Failed replications are set to NA and counted (`cointegration$nfail`) with a warning.
* Monte Carlo size at the 5% level (n = 100, y and x independent random walks, `maxlag = 1`, 200 samples, 199 bootstrap replications; Fov / t / Find, then the combined decision; Monte Carlo standard errors sqrt(r(1 - r)/200) in parentheses). Before = fbardl 1.1.0, after = 1.2.0:
  - Default call (Fourier, k* selected): before (BVZ) 0.375 (0.034) / 0.420 (0.035) / 0.270 (0.031); after (`"fbardl_bvz"`, integer grid, k* and lags selected again) 0.075 (0.019) / 0.135 (0.024) / 0.110 (0.022), decision 0.040 (0.014).
  - `"fbardl_mcnown"` (Fourier, k* selected): before 0.385 (0.034) / 0.435 (0.035) / 0.245 (0.030); after 0.085 (0.020) / 0.130 (0.024) / 0.095 (0.021), decision 0.035 (0.013).
  - `fourier = FALSE`, BVZ: before 0.060 (0.017) / 0.060 (0.017) / 0.125 (0.023); after 0.030 (0.012) / 0.085 (0.020) / 0.085 (0.020), decision 0.020 (0.010).
  - `fourier = FALSE`, MSG: before 0.060 (0.017) / 0.070 (0.018) / 0.085 (0.020); after 0.050 (0.015) / 0.075 (0.019) / 0.080 (0.019), decision 0.035 (0.013).
  - k* = 1 and lags p = 1, q = 1 fixed in advance (`kstar = 1, lags = list(p = 1, q = 1)`), BVZ: before 0.075 (0.019) / 0.095 (0.021) / 0.130 (0.024); after 0.085 (0.020) / 0.110 (0.022) / 0.125 (0.023), decision 0.035 (0.013).
  - The same, MSG: before 0.080 (0.019) / 0.090 (0.020) / 0.105 (0.022); after 0.075 (0.019) / 0.100 (0.021) / 0.105 (0.022), decision 0.025 (0.011).
  - `type = "fardl"`: decision with Fourier terms before 0.605 (0.035), after no decision; without Fourier terms 0.045 (0.015) before and after.
  After the correction Fov, t and Find individually over-reject somewhat (about 0.075 to 0.135 at the 5 percent level with Fourier terms; 0.030 to 0.085 without, for example t 0.135 (0.024) and Find 0.110 (0.022) for the default call, also when k* and the lags are known), while the combined decision does not (0.020 (0.010) to 0.040 (0.014)). The fractional grid was not included in these experiments.
* With `unconditional = TRUE` the bootstrap generating equation for y also omits the contemporaneous differences of the regressors (engine argument `j0 = 1`, the unconditional form of McNown, Sam and Goh, 2018, eq. 12), so the bootstrap model is the estimated model.

## Unchanged results

* For the same inputs (with `kgrid = "fractional"` where the default grid changed), the selected k*, p and q, all coefficient estimates, standard errors, long-run coefficients, diagnostics and the reported Fov, t and Find statistics are identical to fbardl 1.1.0; this was verified bit for bit on 180 combinations of data sets, PSS cases 2 to 4, test types, `fourier` and `unconditional`. Only the bootstrap p-values, critical values and decisions, and the `type = "fardl"` output with Fourier terms, change.

## New arguments

* `kgrid` (`"integer"` or `"fractional"`), `kstar` (Fourier frequency fixed in advance), `lags` (list with `p` and `q` fixed in advance) and `seed` (local seed for the bootstrap; the session's random number generator state is restored).

## Documentation

* The help page states which steps are estimated once on the data and which are redone in every replication, and what is held fixed.
* README: the reference "Bertelli et al., 2020" is corrected to Bertelli, Vacca and Zoia (2022), Economic Modelling 116, 105987, doi:10.1016/j.econmod.2022.105987; the default type and the test types are updated.
* Removed unsupported claims: the Description no longer mentions dynamic multiplier analysis (none is computed; the argument `horizon` is documented as unused), and `unconditional` is no longer attributed to Yilanci, Bozoklu and Gorus (2020), whose paper is not available to us (the simulation model of McNown, Sam and Goh, 2018, eq. 12, has no contemporaneous difference). The duplicated Kripfganz and Schneider reference was removed and Enders and Lee (2012) and Omay (2015) were added.
## Printing

* `summary()`: the significance stars of the diagnostic tests were computed from 1 - p (for example "*" at p = 0.9177); they now mark small p-values (*** p < 0.01, ** p < 0.05, * p < 0.10).
* `summary()` prints the installed package version instead of "v1.0.0".
* `print()` and `summary()` add a note on the numbering of the degenerate cases when one is reported.

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes

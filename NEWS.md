# fbardl 1.1.0

* The PSS case is now applied to the model: cases 4 and 5 include a linear trend, and the overall F test restricts the intercept (case 2) or the trend (case 4) together with the lagged levels. Before, every case estimated the case 3 model.
* `type = "fardl"`: critical values and approximate p-values from the response surfaces of Kripfganz and Schneider (2020), for the sample size, the number of regressors and the number of short-run coefficients. The previous version used asymptotic case 3 values for every case and reported fixed pseudo p-values (0.005, 0.025, 0.075, 0.15, 0.5); the F test on the lagged regressors used an F distribution, which does not apply, and is now reported without critical values.
* Bootstrap types rewritten as a port of the Stata module fbardl 1.3.0: bootstrap data are generated recursively under the null of each test from the restricted equation for y and the equations for Delta x, and the full model is re-estimated on each sample. The previous version did not impose the null (McNown type), imposed only the null of the t test for all three statistics (Bertelli type), and did not regenerate the lagged regressors; its 90% and 95% quantiles were labelled as I(0) and I(1) bounds.
* New argument `unconditional` (no contemporaneous differences of the regressors).
* R-squared is now centred (the model constant is part of the design).
* Breusch-Godfrey, Breusch-Pagan and ARCH diagnostics now follow the standard auxiliary regressions (as in Stata's estat bgodfrey, hettest and archlm); the Breusch-Godfrey test omitted the regressors and the Breusch-Pagan test used a time trend.
* Added tests against Stata reference values.

# fbardl 1.0.3

* Corrected the DOI of Yilanci, Bozoklu and Gorus (2020) to 10.1016/j.scs.2020.102035.
* Authors@R and README author section updated; a former contributor entry was removed.


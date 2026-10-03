# fbardl: Fourier Bootstrap ARDL Cointegration Test

R implementation of the Fourier Bootstrap ARDL (FBARDL) bounds testing approach for cointegration analysis.

## Overview

The `fbardl` package combines:
- **Pesaran, Shin and Smith (2001)** ARDL bounds testing framework
- **Fourier terms** to capture smooth structural breaks (Yilanci et al., 2020)
- **Recursive bootstrap critical values** (McNown et al., 2018; Bertelli et al., 2022), with the Fourier frequency and the lags selected again in every replication (a package choice; or fixed by the user in advance)

The Kripfganz and Schneider (2020) bounds do not account for Fourier terms (in a Monte Carlo with independent random walks the 5% bounds decision rejected in 60.5% of samples), so with Fourier terms inference is by the bootstrap types and the bounds are used for decisions only in models without Fourier terms.

## Installation

```r
install.packages("fbardl")

# Development version
devtools::install_github("muhammedalkhalaf/fbardl")
```

## Usage

```r
library(fbardl)

# Load example data
data(fbardl_data)

# Fourier bootstrap ARDL (Bertelli, Vacca and Zoia scheme, the default);
# k* and the lags are selected again in each bootstrap replication
result <- fbardl(y ~ x1 + x2, data = fbardl_data, reps = 999, seed = 1)
summary(result)

# McNown, Sam and Goh scheme with k* and the lags fixed in advance
result_msg <- fbardl(y ~ x1 + x2, data = fbardl_data, type = "fbardl_mcnown",
                     kstar = 1, lags = list(p = 1, q = c(1, 1)), seed = 1)

# Bounds test (valid without Fourier terms only)
result_bounds <- fbardl(y ~ x1 + x2, data = fbardl_data, type = "fardl",
                        fourier = FALSE)
```

## Test Types

| Type | Description |
|------|-------------|
| `"fbardl_bvz"` | Bootstrap ARDL, separate nulls (Bertelli, Vacca and Zoia, 2022); default |
| `"fbardl_mcnown"` | Bootstrap ARDL, null of the overall F test (McNown, Sam and Goh, 2018) |
| `"fardl"` | Kripfganz and Schneider (2020) bounds; with Fourier terms the bounds are shown for reference only, without a decision |

## Features

- **Automatic lag selection** via AIC or BIC
- **Fourier frequency selection** by minimum SSR over integer frequencies (default) or a fractional grid
- **Three cointegration tests**: F-overall, t-dependent, F-independent
- **Degenerate case detection** (McNown et al., 2018; cases numbered as in Bertelli et al., 2022)
- **Long-run coefficient estimation** via delta method
- **Diagnostic tests**: Jarque-Bera, Breusch-Godfrey, Breusch-Pagan, ARCH

## Output

The function returns an object of class `"fbardl"` containing:
- Model coefficients and standard errors
- Long-run and short-run coefficient estimates
- Cointegration test statistics and p-values
- Diagnostic test results
- Model fit statistics (R-squared, AIC, BIC)

## References

- Bertelli, S., Vacca, G. and Zoia, M. (2022). Bootstrap cointegration tests in ARDL models. *Economic Modelling*, 116, 105987. https://doi.org/10.1016/j.econmod.2022.105987

- Enders, W. and Lee, J. (2012). A unit root test using a Fourier series to approximate smooth breaks. *Oxford Bulletin of Economics and Statistics*, 74(4), 574-599. https://doi.org/10.1111/j.1468-0084.2011.00662.x

- Pesaran, M. H., Shin, Y. and Smith, R. J. (2001). Bounds testing approaches to the analysis of level relationships. *Journal of Applied Econometrics*, 16(3), 289-326. https://doi.org/10.1002/jae.616

- McNown, R., Sam, C. Y. and Goh, S. K. (2018). Bootstrapping the autoregressive distributed lag test for cointegration. *Applied Economics*, 50(13), 1509-1521. https://doi.org/10.1080/00036846.2017.1366643

- Yilanci, V., Bozoklu, S. and Gorus, M. S. (2020). Are BRICS countries pollution havens? Evidence from a bootstrap ARDL bounds testing approach with a Fourier function. *Sustainable Cities and Society*, 55, 102035. https://doi.org/10.1016/j.scs.2020.102035

- Kripfganz, S. and Schneider, D. C. (2020). Response surface regressions for critical value bounds and approximate p-values in equilibrium correction models. *Oxford Bulletin of Economics and Statistics*, 82(6), 1456-1481. https://doi.org/10.1111/obes.12377

## Author

Muhammad Alkhalaf (muhammedalkhalaf@gmail.com)

## License

GPL-3

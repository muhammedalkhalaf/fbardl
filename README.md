# fbardl: Fourier Bootstrap ARDL Cointegration Test

R implementation of the Fourier Bootstrap ARDL (FBARDL) bounds testing approach for cointegration analysis.

## Overview

The `fbardl` package combines:
- **Pesaran, Shin & Smith (2001)** ARDL bounds testing framework
- **Fourier terms** to capture smooth structural breaks (Yilanci et al., 2020)
- **Bootstrap critical values** for robust inference (McNown et al., 2018; Bertelli et al., 2020)

## Installation

```r
# Install from source
install.packages("fbardl_1.0.0.tar.gz", repos = NULL, type = "source")

# Or using devtools
devtools::install_github("muhammedalkhalaf/fbardl")
```

## Usage

```r
library(fbardl)

# Load example data
data(fbardl_data)

# Basic Fourier ARDL test
result <- fbardl(y ~ x1 + x2, data = fbardl_data, type = "fardl")
summary(result)

# Bootstrap ARDL (McNown approach)
result_boot <- fbardl(y ~ x1 + x2, data = fbardl_data,
                      type = "fbardl_mcnown", reps = 999)
summary(result_boot)

# Without Fourier terms
result_nofourier <- fbardl(y ~ x1 + x2, data = fbardl_data, fourier = FALSE)
```

## Test Types

| Type | Description |
|------|-------------|
| `"fardl"` | Fourier ARDL with PSS bounds test (default) |
| `"fbardl_mcnown"` | Bootstrap ARDL (McNown, Sam & Goh, 2018) |
| `"fbardl_bvz"` | Bootstrap ARDL (Bertelli, Vacca & Zoia, 2020) |

## Features

- **Automatic lag selection** via AIC or BIC
- **Fourier frequency selection** by minimum SSR
- **Three cointegration tests**: F-overall, t-dependent, F-independent
- **Degenerate case detection** (McNown et al., 2018)
- **Long-run coefficient estimation** via delta method
- **Diagnostic tests**: Jarque-Bera, Breusch-Godfrey, Breusch-Pagan, ARCH

## Output

The function returns an object of class `"fbardl"` containing:
- Model coefficients and standard errors
- Long-run and short-run coefficient estimates
- Cointegration test statistics and p-values
- Diagnostic test results
- Model fit statistics (R², AIC, BIC)

## References

- Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches to the analysis of level relationships. *Journal of Applied Econometrics*, 16(3), 289-326. https://doi.org/10.1002/jae.616

- McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the autoregressive distributed lag test for cointegration. *Applied Economics*, 50(13), 1509-1521. https://doi.org/10.1080/00036846.2017.1366643

- Yilanci, V., Bozoklu, S., & Gorus, M. S. (2020). Are BRICS countries pollution havens? Evidence from a bootstrap ARDL bounds testing approach with a Fourier function. *Sustainable Cities and Society*, 55, 102035. https://doi.org/10.1016/j.scs.2020.102035

- Kripfganz, S., & Schneider, D. C. (2020). Response surface regressions for critical value bounds and approximate p-values in equilibrium correction models. *Oxford Bulletin of Economics and Statistics*, 82(6), 1456-1481. https://doi.org/10.1111/obes.12377

## Author

Dr. Merwan Roudane (merwanroudane920@gmail.com)

## License

GPL-3

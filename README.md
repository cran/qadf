# qadf

**Quantile ADF Unit Root Test** for R

## Overview

`qadf` implements the Quantile Autoregressive Distributed Lag (QADF) unit
root test proposed by Koenker and Xiao (2004, JASA). The test uses quantile
regression to examine unit root behaviour across the conditional distribution
of a time series, providing a richer characterisation of persistence than
classical ADF tests.

## Installation

```r
# From CRAN (once published):
install.packages("qadf")
```

## Usage

```r
library(qadf)

set.seed(42)
y <- cumsum(rnorm(120))

# Test at median (tau = 0.5) with constant model
result <- qadf(y, tau = 0.5, model = "c", max_lags = 8, ic = "aic")
print(result)

# Test at lower tail (tau = 0.25) with trend model
result_q25 <- qadf(y, tau = 0.25, model = "ct", ic = "bic")
print(result_q25)
```

## References

Koenker, R. and Xiao, Z. (2004). Unit Root Quantile Autoregression Inference.
*Journal of the American Statistical Association*, 99(465), 775–787.
<https://doi.org/10.1198/016214504000001114>

Hansen, B. E. (1995). Rethinking the Univariate Approach to Unit Root Tests.
*Econometric Theory*, 11(5), 1148–1171.
<https://doi.org/10.1017/S0266466600009713>

## License

GPL-3

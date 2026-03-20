\name{qadf}
\alias{qadf}
\title{Quantile ADF Unit Root Test}
\usage{
qadf(x, tau = 0.5, model = "c", max_lags = 8, ic = "aic")
}
\arguments{
\item{x}{A numeric vector or univariate time series object.}

\item{tau}{A numeric scalar specifying the quantile at which to estimate
the model. Must satisfy \code{0 < tau < 1}. Default is \code{0.5}.}

\item{model}{A character string specifying the deterministic component.
\code{"c"} (default) includes a constant; \code{"ct"} includes a
constant and a linear trend.}

\item{max_lags}{A non-negative integer specifying the maximum number of
augmentation lags to consider. Default is \code{8}.}

\item{ic}{A character string for the information criterion used to select
the optimal lag length. One of \code{"aic"} (default), \code{"bic"}, or
\code{"tstat"} (sequential t-test at the 10\% level).}
}
\value{
An object of class \code{"qadf"} with components:
\describe{
  \item{statistic}{The QADF t-statistic \eqn{t_n(\tau)}.}
  \item{coef_stat}{The \eqn{U_n(\tau) = n(\hat\rho(\tau) - 1)} statistic.}
  \item{rho_tau}{Quantile autoregressive coefficient \eqn{\hat\rho(\tau)}.}
  \item{rho_ols}{OLS autoregressive coefficient.}
  \item{alpha_tau}{Quantile intercept \eqn{\hat\alpha_0(\tau)}.}
  \item{delta2}{Nuisance parameter \eqn{\hat\delta^2}.}
  \item{half_life}{Half-life implied by \eqn{\hat\rho(\tau)}, in periods.}
  \item{opt_lags}{Selected lag order.}
  \item{nobs}{Number of observations used.}
  \item{critical_values}{Named numeric vector of critical values at 1\%,
    5\%, and 10\% from Hansen (1995).}
  \item{tau}{The quantile used.}
  \item{model}{The deterministic model used.}
  \item{ic}{The information criterion used.}
  \item{varname}{The name of the input series.}
}
}
\description{
Implements the Quantile Autoregressive Distributed Lag (QADF) unit root
test of Koenker and Xiao (2004). The test examines unit root behaviour
across quantiles of the conditional distribution of a time series using
quantile regression.
}
\details{
The QADF test estimates the autoregressive parameter \eqn{\hat\rho(\tau)}
at quantile \eqn{\tau} via quantile regression on the ADF regression
equation. The t-statistic \eqn{t_n(\tau) = (\hat\rho(\tau) - 1) / se}
tests \eqn{H_0: \rho(\tau) = 1} (unit root) against
\eqn{H_1: \rho(\tau) < 1} (stationarity).

Critical values are from Table 1 of Hansen (1995), interpolated
linearly for quantiles between tabulated values. The model
\code{"c"} corresponds to a demeaned ADF regression; \code{"ct"}
adds a linear time trend.
}
\examples{
set.seed(42)
y <- cumsum(rnorm(100))
result <- qadf(y, tau = 0.5, model = "c", max_lags = 4)
print(result)
}
\references{
Koenker, R. and Xiao, Z. (2004). Unit Root Quantile Autoregression
Inference. \emph{Journal of the American Statistical Association},
99(465), 775--787. \doi{10.1198/016214504000001114}

Hansen, B. E. (1995). Rethinking the Univariate Approach to Unit Root
Tests: How to Use Covariates to Increase Power. \emph{Econometric Theory},
11(5), 1148--1171. \doi{10.1017/S0266466600009713}
}

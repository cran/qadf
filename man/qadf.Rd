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
\code{"tstat"} (largest lag whose last coefficient is significant at the 5\% level).}
}
\value{
An object of class \code{"qadf"} with components:
\describe{
  \item{statistic}{The QADF t-statistic \eqn{t_n(\tau)}.}
  \item{coef_stat}{The \eqn{U_n(\tau) = n(\hat\rho(\tau) - 1)} statistic.}
  \item{rho_tau}{Quantile autoregressive coefficient \eqn{\hat\rho(\tau)}.}
  \item{rho_ols}{OLS autoregressive coefficient.}
  \item{alpha_tau}{Quantile intercept \eqn{\hat\alpha_0(\tau)}.}
  \item{delta2}{Nuisance parameter \eqn{\hat\delta^2}, the squared
    correlation between \eqn{\Delta y_t} and \eqn{\psi_\tau(\hat u_t)};
    it indexes the critical values.}
  \item{half_life}{Half-life implied by \eqn{\hat\rho(\tau)}, in periods.}
  \item{opt_lags}{Selected lag order.}
  \item{nobs}{Number of observations used.}
  \item{critical_values}{Named numeric vector of critical values at 1\%,
    5\%, and 10\% from Hansen (1995), interpolated in \eqn{\hat\delta^2}.}
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
The quantile autoregression is estimated in levels,
\deqn{Q_{y_t}(\tau) = \alpha_0(\tau) + \rho(\tau) y_{t-1} +
  \sum_{j=1}^{p} \alpha_j(\tau) \Delta y_{t-j},}
with a linear trend added for \code{model = "ct"}. The lag order \eqn{p}
is chosen from the ADF regression on a common sample. The statistic is
equation (9) of Koenker and Xiao (2004),
\deqn{t_n(\tau) = \frac{\hat f(F^{-1}(\tau))}{\sqrt{\tau(1-\tau)}}
  (Y_{-1}' P_X Y_{-1})^{1/2} (\hat\rho(\tau) - 1),}
where \eqn{P_X} projects off the other regressors and the density is
estimated by the difference quotient of the fitted conditional quantile at
\eqn{\tau \pm h} with the Hall-Sheather bandwidth \eqn{h}. It tests
\eqn{H_0: \rho(\tau) = 1} against \eqn{H_1: \rho(\tau) < 1}.

The limiting distribution depends on \eqn{\delta^2}; the critical values
are those of Hansen (1995) for the covariate augmented Dickey-Fuller test,
linearly interpolated on the grid \eqn{\delta^2 = 0.1, \ldots, 1}.
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
11(5), 1148--1171. \doi{10.1017/S0266466600009993}
}

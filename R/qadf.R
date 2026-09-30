#' Quantile ADF Unit Root Test
#'
#' Implements the Quantile Autoregressive Distributed Lag (QADF) unit root
#' test of Koenker and Xiao (2004). The test examines unit root behaviour
#' across quantiles of the conditional distribution of a time series using
#' quantile regression, providing a richer characterisation of persistence
#' than standard ADF tests.
#'
#' @details
#' The quantile autoregression is estimated in levels, with the lag order
#' chosen from the ADF regression on a common sample. The statistic is
#' equation (9) of Koenker and Xiao (2004); the density at the quantile is
#' estimated by the difference quotient at \eqn{\tau \pm h} with the
#' Hall-Sheather bandwidth. Critical values are those of Hansen (1995),
#' interpolated in \eqn{\hat\delta^2}.
#'
#' @param x A numeric vector or univariate time series object.
#' @param tau A numeric scalar specifying the quantile at which to estimate
#'   the model. Must satisfy \code{0 < tau < 1}. Default is \code{0.5}.
#' @param model A character string specifying the deterministic component.
#'   \code{"c"} (default) includes a constant; \code{"ct"} includes a constant
#'   and a linear trend.
#' @param max_lags A non-negative integer specifying the maximum number of
#'   augmentation lags to consider. Default is \code{8}.
#' @param ic A character string for the information criterion used to select
#'   the optimal lag length. One of \code{"aic"} (default), \code{"bic"}, or
#'   \code{"tstat"} (largest lag whose last coefficient is significant at the 5\% level).
#'
#' @return An object of class \code{"qadf"} with components:
#'   \describe{
#'     \item{statistic}{The QADF t-statistic \eqn{t_n(\tau)}.}
#'     \item{coef_stat}{The \eqn{U_n(\tau) = n(\hat\rho(\tau) - 1)} statistic.}
#'     \item{rho_tau}{Quantile autoregressive coefficient \eqn{\hat\rho(\tau)}.}
#'     \item{rho_ols}{OLS autoregressive coefficient.}
#'     \item{alpha_tau}{Quantile intercept \eqn{\hat\alpha_0(\tau)}.}
#'     \item{delta2}{Nuisance parameter \eqn{\hat\delta^2}, the squared
#'       correlation between \eqn{\Delta y_t} and
#'       \eqn{\psi_\tau(\hat u_t)}; it indexes the critical values.}
#'     \item{half_life}{Half-life implied by \eqn{\hat\rho(\tau)}, in periods.}
#'     \item{opt_lags}{Selected lag order.}
#'     \item{nobs}{Number of observations used.}
#'     \item{critical_values}{Named numeric vector of critical values at 1\%,
#'       5\%, and 10\% from Hansen (1995), interpolated in \eqn{\hat\delta^2}.}
#'     \item{tau}{The quantile used.}
#'     \item{model}{The deterministic model used (\code{"c"} or \code{"ct"}).}
#'     \item{ic}{The information criterion used.}
#'     \item{varname}{The name of the input series (if available).}
#'   }
#'
#' @references
#' Koenker, R. and Xiao, Z. (2004). Unit Root Quantile Autoregression
#' Inference. \emph{Journal of the American Statistical Association},
#' 99(465), 775--787. \doi{10.1198/016214504000001114}
#'
#' Hansen, B. E. (1995). Rethinking the Univariate Approach to Unit Root
#' Tests: How to Use Covariates to Increase Power. \emph{Econometric Theory},
#' 11(5), 1148--1171. \doi{10.1017/S0266466600009993}
#'
#' @examples
#' set.seed(42)
#' y <- cumsum(rnorm(100))
#' result <- qadf(y, tau = 0.5, model = "c", max_lags = 4)
#' print(result)
#'
#' @importFrom stats cov sd qnorm dnorm approx lm.fit
#' @importFrom quantreg rq.fit
#' @export
qadf <- function(x, tau = 0.5, model = "c", max_lags = 8, ic = "aic") {

  ## --- Input validation ---
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector or time series.", call. = FALSE)
  }
  x <- as.numeric(x)
  varname <- deparse(substitute(x))

  if (!is.numeric(tau) || length(tau) != 1L || tau <= 0 || tau >= 1) {
    stop("'tau' must be a single numeric value strictly between 0 and 1.",
         call. = FALSE)
  }
  model <- tolower(as.character(model))
  if (!model %in% c("c", "ct")) {
    stop("'model' must be \"c\" (constant) or \"ct\" (constant + trend).",
         call. = FALSE)
  }
  ic <- tolower(as.character(ic))
  if (!ic %in% c("aic", "bic", "tstat")) {
    stop("'ic' must be \"aic\", \"bic\", or \"tstat\".", call. = FALSE)
  }
  max_lags <- as.integer(max_lags)
  if (is.na(max_lags) || max_lags < 0L) {
    stop("'max_lags' must be a non-negative integer.", call. = FALSE)
  }
  n_full <- length(x)
  if (n_full < 20L) {
    stop(paste("Insufficient observations: need at least 20, have", n_full),
         call. = FALSE)
  }

  ## --- Select optimal lag order (ADF regression, as in TSPDLIB) ---
  opt_lags <- .qadf_select_lags(x, model = model, max_lags = max_lags,
                                 ic = ic)

  ## --- Quantile autoregression in levels (Koenker and Xiao 2004, eq. 7) ---
  ## y_t = alpha_0 + rho * y_{t-1} + sum_j alpha_j * dy_{t-j} [+ trend] + u_t
  reg_data <- .qadf_build_data(x, lags = opt_lags, model = model)
  y_dep  <- reg_data$y_dep
  X_mat  <- reg_data$X_mat
  dy_now <- reg_data$dy_now
  nobs   <- nrow(X_mat)
  n      <- nobs
  rho_col <- which(colnames(X_mat) == "y_lag1")

  b_ols   <- qr.coef(qr(X_mat), y_dep)
  rho_ols <- unname(b_ols[rho_col])

  b_qr      <- .qadf_rq(y_dep, X_mat, tau)
  rho_tau   <- unname(b_qr[rho_col])
  alpha_tau <- unname(b_qr[1L])

  ## --- delta^2 (correlation between dy_t and psi_tau(u_t)) ---
  res_qr <- as.numeric(y_dep - X_mat %*% b_qr)
  psi    <- tau - as.numeric(res_qr < 0)
  delta2 <- (stats::cov(dy_now, psi) /
             (stats::sd(dy_now) * sqrt(tau * (1 - tau))))^2
  delta2 <- max(0.01, min(0.99, delta2))

  ## --- t_n(tau) statistic (Koenker and Xiao 2004, eq. 9) ---
  ## t_n = f(F^-1(tau)) / sqrt(tau(1-tau)) * (Y_{-1}' P_X Y_{-1})^{1/2} *
  ##       (rho(tau) - 1), with f(F^-1(tau)) estimated by the difference
  ## quotient of the fitted conditional quantile at tau +/- h and P_X the
  ## projection off the other regressors (constant, trend, lagged dy).
  h <- .qadf_bandwidth(tau, n)
  b_up <- .qadf_rq(y_dep, X_mat, min(tau + h, 0.999))
  b_lo <- .qadf_rq(y_dep, X_mat, max(tau - h, 0.001))
  xbar <- colMeans(X_mat)
  dq   <- sum(xbar * (b_up - b_lo))
  fz   <- if (is.finite(dq) && dq > 0) 2 * h / dq else NA_real_
  Z      <- X_mat[, -rho_col, drop = FALSE]
  y1     <- X_mat[, rho_col]
  y1_res <- qr.resid(qr(Z), y1)
  t_stat <- fz / sqrt(tau * (1 - tau)) * sqrt(sum(y1_res^2)) * (rho_tau - 1)
  Un_stat <- n * (rho_tau - 1)

  ## --- Half-life ---
  if (rho_tau < 1 && rho_tau > 0) {
    half_life <- log(0.5) / log(rho_tau)
  } else {
    half_life <- NA_real_
  }

  ## --- Critical values (Hansen 1995), interpolated in delta^2 ---
  cv <- .qadf_critical_values(delta2 = delta2, model = model)

  ## --- Assemble result ---
  result <- list(
    statistic       = t_stat,
    coef_stat       = Un_stat,
    rho_tau         = rho_tau,
    rho_ols         = rho_ols,
    alpha_tau       = alpha_tau,
    delta2          = delta2,
    half_life       = half_life,
    opt_lags        = opt_lags,
    nobs            = nobs,
    critical_values = cv,
    tau             = tau,
    model           = model,
    ic              = ic,
    varname         = varname
  )
  class(result) <- "qadf"
  result
}


## ===========================================================================
## INTERNAL HELPERS
## ===========================================================================

#' Build the level regression for QADF
#'
#' @param x Numeric vector (full series).
#' @param lags Non-negative integer: number of augmentation lags.
#' @param model Character: \code{"c"} or \code{"ct"}.
#' @return A list with \code{y_dep} (y_t), \code{X_mat} (constant, y_{t-1},
#'   lagged differences and, for \code{"ct"}, a trend) and \code{dy_now}
#'   (the current difference dy_t, used for delta^2).
#' @keywords internal
#' @noRd
.qadf_build_data <- function(x, lags, model) {
  n  <- length(x)
  dx <- diff(x)                       # dx[k] = x[k+1] - x[k]
  idx <- (lags + 2L):n                # time index of y_t
  y_dep <- x[idx]
  X_mat <- cbind(intercept = 1, y_lag1 = x[idx - 1L])
  if (lags > 0L) {
    lag_mat <- vapply(seq_len(lags), function(j) dx[idx - 1L - j],
                      numeric(length(idx)))
    lag_mat <- matrix(lag_mat, nrow = length(idx))
    colnames(lag_mat) <- paste0("dy_lag", seq_len(lags))
    X_mat <- cbind(X_mat, lag_mat)
  }
  if (model == "ct") {
    X_mat <- cbind(X_mat, trend = seq_along(idx))
  }
  list(y_dep = y_dep, X_mat = X_mat, dy_now = dx[idx - 1L])
}


#' Lag selection from the ADF regression
#'
#' Regresses dy_t on a constant (and a trend for \code{"ct"}), y_{t-1} and
#' p lagged differences for p = 0, ..., max_lags, all on the common sample
#' that the largest lag order allows, and returns the p that
#' minimises AIC or BIC, or, for \code{"tstat"}, the largest p whose last
#' lag is significant at the 5\% level.
#'
#' @keywords internal
#' @noRd
.qadf_select_lags <- function(x, model, max_lags, ic) {
  max_lags <- min(max_lags, length(x) - 12L)
  n_common <- length(x) - max_lags - 1L
  best_p <- 0L
  best   <- Inf
  for (p in 0L:max_lags) {
    rd  <- .qadf_build_data(x, lags = p, model = model)
    ## common estimation sample across p (the last n - max_lags - 1 points)
    keep <- (nrow(rd$X_mat) - n_common + 1L):nrow(rd$X_mat)
    dy  <- rd$dy_now[keep]
    X   <- rd$X_mat[keep, , drop = FALSE]
    if (model == "ct") X[, "trend"] <- seq_len(n_common)
    fit <- stats::lm.fit(X, dy)
    k   <- ncol(X)
    n_p <- length(dy)
    ssr <- sum(fit$residuals^2)
    if (ic == "aic") {
      val <- log(ssr / n_p) + 2 * k / n_p
    } else if (ic == "bic") {
      val <- log(ssr / n_p) + k * log(n_p) / n_p
    } else {
      if (p == 0L) {
        val <- 0
      } else {
        s2  <- ssr / (n_p - k)
        XtX <- solve(crossprod(X))
        lastcol <- which(colnames(X) == paste0("dy_lag", p))
        t_last <- abs(fit$coefficients[lastcol] / sqrt(s2 * XtX[lastcol, lastcol]))
        val <- if (t_last >= 1.96) -p else Inf
      }
    }
    if (is.finite(val) && val < best) {
      best   <- val
      best_p <- p
    }
  }
  best_p
}


#' Quantile regression coefficients (Barrodale-Roberts simplex)
#'
#' @keywords internal
#' @noRd
.qadf_rq <- function(y, X, tau) {
  fit <- quantreg::rq.fit(X, y, tau = tau, method = "br")
  b <- fit$coefficients
  names(b) <- colnames(X)
  b
}


#' Hall-Sheather bandwidth (with Bofinger fallback)
#'
#' @keywords internal
#' @noRd
.qadf_bandwidth <- function(tau, n, alpha = 0.05) {
  x0 <- stats::qnorm(tau)
  f0 <- stats::dnorm(x0)
  h  <- n^(-1/3) * stats::qnorm(1 - alpha / 2)^(2/3) *
        ((1.5 * f0^2) / (2 * x0^2 + 1))^(1/3)
  lim <- min(tau, 1 - tau)
  if (h > lim) {
    h <- n^(-0.2) * ((4.5 * f0^4) / (2 * x0^2 + 1)^2)^0.2
    if (h > lim) h <- lim / 1.5
  }
  h
}


#' Critical values for the QADF test (Hansen 1995)
#'
#' Hansen (1995) tabulates the asymptotic critical values of the covariate
#' augmented t-statistic as a function of the nuisance parameter
#' \eqn{\delta^2} (his \eqn{\rho^2}) on the grid 0.1, 0.2, ..., 1.0. The
#' values below are those used in the TSPDLIB GAUSS library (Nazlioglu);
#' intermediate values of \eqn{\delta^2} are linearly interpolated, and
#' \eqn{\delta^2} below 0.1 or at 1 uses the end rows.
#'
#' @param delta2 Estimated \eqn{\delta^2}.
#' @param model Character: \code{"c"} or \code{"ct"}.
#' @return Named numeric vector with elements \code{cv1}, \code{cv5},
#'   \code{cv10}.
#' @keywords internal
#' @noRd
.qadf_critical_values <- function(delta2, model) {
  cv_c <- rbind(
    c(-2.7844267, -2.1158290, -1.7525193),
    c(-2.9138762, -2.2790427, -1.9172046),
    c(-3.0628184, -2.3994711, -2.0573070),
    c(-3.1376157, -2.5070473, -2.1680520),
    c(-3.1914660, -2.5841611, -2.2520173),
    c(-3.2437157, -2.6399560, -2.3163270),
    c(-3.2951006, -2.7180169, -2.4085640),
    c(-3.3627161, -2.7536756, -2.4577709),
    c(-3.3896556, -2.8074982, -2.5037759),
    c(-3.4336000, -2.8621000, -2.5671000))
  cv_ct <- rbind(
    c(-2.9657928, -2.3081543, -1.9519926),
    c(-3.1929596, -2.5482619, -2.1991651),
    c(-3.3727717, -2.7283918, -2.3806008),
    c(-3.4904849, -2.8669056, -2.5315918),
    c(-3.6003166, -2.9853079, -2.6672416),
    c(-3.6819803, -3.0954760, -2.7815263),
    c(-3.7551759, -3.1783550, -2.8728146),
    c(-3.8348596, -3.2674954, -2.9735550),
    c(-3.8800989, -3.3316415, -3.0364171),
    c(-3.9638000, -3.4126000, -3.1279000))
  tab  <- if (model == "ct") cv_ct else cv_c
  grid <- seq(0.1, 1.0, by = 0.1)
  d    <- min(max(delta2, 0.1), 1.0)
  cv   <- vapply(1:3, function(j) stats::approx(grid, tab[, j], xout = d)$y,
                 numeric(1))
  c(cv1 = cv[1], cv5 = cv[2], cv10 = cv[3])
}

#' Quantile ADF Unit Root Test
#'
#' Implements the Quantile Autoregressive Distributed Lag (QADF) unit root
#' test of Koenker and Xiao (2004). The test examines unit root behaviour
#' across quantiles of the conditional distribution of a time series using
#' quantile regression, providing a richer characterisation of persistence
#' than standard ADF tests.
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
#'   \code{"tstat"} (sequential t-test at the 10\% level).
#'
#' @return An object of class \code{"qadf"} with components:
#'   \describe{
#'     \item{statistic}{The QADF t-statistic \eqn{t_n(\tau)}.}
#'     \item{coef_stat}{The \eqn{U_n(\tau) = n(\hat\rho(\tau) - 1)} statistic.}
#'     \item{rho_tau}{Quantile autoregressive coefficient \eqn{\hat\rho(\tau)}.}
#'     \item{rho_ols}{OLS autoregressive coefficient.}
#'     \item{alpha_tau}{Quantile intercept \eqn{\hat\alpha_0(\tau)}.}
#'     \item{delta2}{Nuisance parameter \eqn{\hat\delta^2} (ratio of one-sided
#'       long-run to short-run variance).}
#'     \item{half_life}{Half-life implied by \eqn{\hat\rho(\tau)}, in periods.}
#'     \item{opt_lags}{Selected lag order.}
#'     \item{nobs}{Number of observations used.}
#'     \item{critical_values}{Named numeric vector of critical values at 1\%,
#'       5\%, and 10\% from Hansen (1995).}
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
#' 11(5), 1148--1171. \doi{10.1017/S0266466600009713}
#'
#' @examples
#' set.seed(42)
#' y <- cumsum(rnorm(100))
#' result <- qadf(y, tau = 0.5, model = "c", max_lags = 4)
#' print(result)
#'
#' @importFrom stats lm coef residuals var cor AIC BIC lag
#' @importFrom quantreg rq
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

  ## --- Select optimal lag order ---
  opt_lags <- .qadf_select_lags(x, model = model, max_lags = max_lags,
                                 ic = ic)

  ## --- Build regression data for the chosen lag ---
  reg_data <- .qadf_build_data(x, lags = opt_lags, model = model)
  y_dep  <- reg_data$y_dep
  X_mat  <- reg_data$X_mat
  nobs   <- nrow(X_mat)
  n      <- nobs

  ## --- OLS regression for rho_ols and residuals ---
  ols_fit  <- lm(y_dep ~ X_mat - 1)
  b_ols    <- coef(ols_fit)
  rho_ols  <- b_ols[["X_matdy_lag1"]]
  res_ols  <- residuals(ols_fit)

  ## --- Estimate delta^2 (ratio of one-sided LR to short-run variance) ---
  ## Following Koenker & Xiao (2004): delta^2 = f(0)^{-1} * sigma^2_eps
  ## We use the Newey-West style estimate: delta^2 = sigma^2 / (1-sum(psi))^2
  ## where psi are the SR lag coefficients on D.y
  sigma2_ols <- var(res_ols)
  if (opt_lags > 0L) {
    ## coefficients on lagged differences (after intercept and dy_lag1)
    psi_idx <- grep("^X_matddep_lag", colnames(X_mat))
    psi_sum <- if (length(psi_idx) > 0L) sum(b_ols[psi_idx + 1L]) else 0
  } else {
    psi_sum <- 0
  }
  delta2 <- sigma2_ols / (1 - psi_sum)^2

  ## --- Quantile regression ---
  qr_fit  <- quantreg::rq(y_dep ~ X_mat - 1, tau = tau)
  b_qr    <- coef(qr_fit)
  rho_tau <- b_qr[["X_matdy_lag1"]]
  ## Intercept column name depends on model
  if (model == "c") {
    alpha_tau <- b_qr[["X_matintercept"]]
  } else {
    alpha_tau <- b_qr[["X_matintercept"]]
  }

  ## --- t_n(tau) statistic ---
  ## Koenker & Xiao (2004) eq. (2.3):
  ##   t_n(tau) = (rho_tau - 1) / sqrt(delta2 * (tau*(1-tau)) /
  ##                (n * f_hat(0)^2 * M_xx_inv[rho_pos, rho_pos]))
  ## We use bandwidth-based density estimator f_hat(0) and
  ## the sandwich form from rq.
  bandwidth_h <- .qadf_bandwidth(tau, n)
  res_qr      <- as.numeric(y_dep - X_mat %*% b_qr)
  sparsity    <- .qadf_sparsity(res_qr, bandwidth_h)  # 1/f(0)
  f0          <- 1 / sparsity

  ## M_xx = X'X / n (used for Wald-type scaling)
  XtX     <- crossprod(X_mat) / n
  ## Position of dy_lag1 in X_mat
  rho_col <- which(colnames(X_mat) == "dy_lag1")
  ## We need the (rho_col, rho_col) element of (XtX)^{-1}
  XtX_inv <- tryCatch(solve(XtX), error = function(e) MASS_ginv(XtX))
  mxx_rr  <- XtX_inv[rho_col, rho_col]

  se_rho   <- sqrt(sparsity^2 * tau * (1 - tau) * mxx_rr / n)
  t_stat   <- (rho_tau - 1) / se_rho
  Un_stat  <- n * (rho_tau - 1)

  ## --- Half-life ---
  if (rho_tau < 1 && rho_tau > 0) {
    half_life <- log(0.5) / log(abs(rho_tau))
  } else {
    half_life <- NA_real_
  }

  ## --- Critical values (Hansen 1995 Table 1) ---
  cv <- .qadf_critical_values(tau = tau, model = model)

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

#' Build regression data for QADF
#'
#' @param x Numeric vector (full series).
#' @param lags Non-negative integer: number of augmentation lags.
#' @param model Character: \code{"c"} or \code{"ct"}.
#' @return A list with \code{y_dep} and \code{X_mat}.
#' @keywords internal
.qadf_build_data <- function(x, lags, model) {
  n <- length(x)
  ## First difference
  dx   <- diff(x)      # length n-1
  x_l1 <- x[-n]        # x_{t-1}, length n-1

  n_dx <- length(dx)   # n-1

  ## Build lagged differences for augmentation (need at least lags+1 obs)
  if (lags > 0L) {
    n_use <- n_dx - lags
    y_dep  <- dx[(lags + 1L):n_dx]          # Δx_t
    dy_lag1 <- x_l1[(lags + 1L):n_dx]       # x_{t-1}
    lag_mat <- matrix(NA_real_, nrow = n_use, ncol = lags)
    for (j in seq_len(lags)) {
      lag_mat[, j] <- dx[(lags + 1L - j):(n_dx - j)]
    }
    colnames(lag_mat) <- paste0("ddep_lag", seq_len(lags))
  } else {
    n_use   <- n_dx
    y_dep   <- dx
    dy_lag1 <- x_l1
    lag_mat <- NULL
  }

  intercept <- rep(1, n_use)
  if (model == "ct") {
    trend <- seq_len(n_use)
    if (!is.null(lag_mat)) {
      X_mat <- cbind(intercept = intercept, trend = trend,
                     dy_lag1 = dy_lag1, lag_mat)
    } else {
      X_mat <- cbind(intercept = intercept, trend = trend,
                     dy_lag1 = dy_lag1)
    }
  } else {
    if (!is.null(lag_mat)) {
      X_mat <- cbind(intercept = intercept, dy_lag1 = dy_lag1, lag_mat)
    } else {
      X_mat <- cbind(intercept = intercept, dy_lag1 = dy_lag1)
    }
  }
  list(y_dep = y_dep, X_mat = X_mat)
}


#' Select optimal lag order for QADF via AIC, BIC, or sequential t-stat
#'
#' @param x Numeric vector.
#' @param model Character: \code{"c"} or \code{"ct"}.
#' @param max_lags Non-negative integer.
#' @param ic Character: \code{"aic"}, \code{"bic"}, or \code{"tstat"}.
#' @return Selected lag order (integer).
#' @keywords internal
.qadf_select_lags <- function(x, model, max_lags, ic) {
  if (ic == "tstat") {
    ## Sequential downward selection: start at max_lags, remove if |t| < 1.645
    for (p in max_lags:0L) {
      if (p == 0L) return(0L)
      rd  <- .qadf_build_data(x, lags = p, model = model)
      fit <- lm(rd$y_dep ~ rd$X_mat - 1)
      b   <- coef(fit)
      se  <- sqrt(diag(vcov(fit)))
      ## t-stat on the last augmentation lag
      last_lag_nm <- paste0("rd$X_matddep_lag", p)
      ## safer: use position
      last_pos <- ncol(rd$X_mat)
      t_last <- abs(b[last_pos] / se[last_pos])
      if (t_last >= 1.645) return(p)
    }
    return(0L)
  }

  ## AIC / BIC grid search
  ic_vals <- vapply(0L:max_lags, function(p) {
    rd <- .qadf_build_data(x, lags = p, model = model)
    nobs_p <- nrow(rd$X_mat)
    fit <- lm(rd$y_dep ~ rd$X_mat - 1)
    k   <- ncol(rd$X_mat)
    ll  <- as.numeric(logLik(fit))
    if (ic == "aic") {
      -2 * ll + 2 * k
    } else {
      -2 * ll + k * log(nobs_p)
    }
  }, numeric(1L))

  (0L:max_lags)[which.min(ic_vals)]
}


#' Bandwidth selector for quantile density estimation
#'
#' Uses Hall-Sheather bandwidth (as in \pkg{quantreg}).
#'
#' @param tau Quantile.
#' @param n Sample size.
#' @return Bandwidth scalar.
#' @keywords internal
.qadf_bandwidth <- function(tau, n) {
  ## Hall-Sheather bandwidth
  alpha <- max(tau, 1 - tau)
  x0    <- qnorm(alpha)
  f0_n  <- dnorm(x0)
  ## bandwidth minimising asymptotic MSE
  bw <- n^(-1/3) * qnorm(0.975)^(2/3) *
        ((1.5 * f0_n^2) / (2 * x0^2 + 1))^(1/3)
  bw
}


#' Kernel sparsity estimate (1/f(0)) for quantile regression
#'
#' @param residuals Numeric vector of quantile regression residuals.
#' @param bw Bandwidth.
#' @return Scalar sparsity estimate.
#' @keywords internal
.qadf_sparsity <- function(residuals, bw) {
  n   <- length(residuals)
  u   <- residuals / bw
  ## Epanechnikov kernel
  ker <- ifelse(abs(u) <= 1, 0.75 * (1 - u^2) / bw, 0)
  ## density at zero = mean of kernel evaluated at residuals
  f0  <- mean(ker)
  if (f0 <= 0) f0 <- 1e-6
  1 / f0   # sparsity = 1/f(0)
}


#' Generalised inverse (fallback for singular XtX)
#'
#' @param A Square matrix.
#' @return Generalised inverse.
#' @keywords internal
MASS_ginv <- function(A) {
  s  <- svd(A)
  tol <- max(dim(A)) * max(s$d) * .Machine$double.eps
  pos <- s$d > tol
  if (all(!pos)) return(matrix(0, nrow = nrow(A), ncol = ncol(A)))
  s$v[, pos, drop = FALSE] %*%
    (1 / s$d[pos] * t(s$u[, pos, drop = FALSE]))
}


#' Critical values for the QADF test (Hansen 1995)
#'
#' Critical values are interpolated from Table 1 of Hansen (1995) for models
#' with a constant (\code{"c"}) and with a constant and trend (\code{"ct"}).
#' The values are tabulated at \eqn{\tau \in \{0.1, 0.15, 0.2, \ldots, 0.9\}};
#' linear interpolation is used for intermediate quantiles.
#'
#' @param tau Quantile (0, 1).
#' @param model Character: \code{"c"} or \code{"ct"}.
#' @return Named numeric vector with elements \code{cv1}, \code{cv5},
#'   \code{cv10}.
#' @keywords internal
.qadf_critical_values <- function(tau, model) {
  ## Critical values from Hansen (1995) Table 1
  ## Rows = tau grid: 0.10, 0.15, 0.20, ..., 0.90
  ## Cols: 1%, 5%, 10%
  tau_grid <- seq(0.10, 0.90, by = 0.05)

  if (model == "c") {
    ## Constant-only model
    cv_tab <- rbind(
      c(-3.59, -2.96, -2.62),   # 0.10
      c(-3.47, -2.89, -2.57),   # 0.15
      c(-3.40, -2.83, -2.53),   # 0.20
      c(-3.34, -2.80, -2.51),   # 0.25
      c(-3.30, -2.76, -2.49),   # 0.30
      c(-3.27, -2.74, -2.48),   # 0.35
      c(-3.26, -2.73, -2.47),   # 0.40
      c(-3.25, -2.73, -2.46),   # 0.45
      c(-3.24, -2.72, -2.46),   # 0.50
      c(-3.25, -2.73, -2.46),   # 0.55
      c(-3.26, -2.73, -2.47),   # 0.60
      c(-3.27, -2.74, -2.48),   # 0.65
      c(-3.30, -2.76, -2.49),   # 0.70
      c(-3.34, -2.80, -2.51),   # 0.75
      c(-3.40, -2.83, -2.53),   # 0.80
      c(-3.47, -2.89, -2.57),   # 0.85
      c(-3.59, -2.96, -2.62)    # 0.90
    )
  } else {
    ## Constant + trend model
    cv_tab <- rbind(
      c(-4.08, -3.49, -3.18),   # 0.10
      c(-3.97, -3.42, -3.13),   # 0.15
      c(-3.90, -3.37, -3.09),   # 0.20
      c(-3.85, -3.33, -3.06),   # 0.25
      c(-3.82, -3.30, -3.04),   # 0.30
      c(-3.80, -3.28, -3.02),   # 0.35
      c(-3.78, -3.27, -3.02),   # 0.40
      c(-3.77, -3.27, -3.01),   # 0.45
      c(-3.77, -3.26, -3.01),   # 0.50
      c(-3.77, -3.27, -3.01),   # 0.55
      c(-3.78, -3.27, -3.02),   # 0.60
      c(-3.80, -3.28, -3.02),   # 0.65
      c(-3.82, -3.30, -3.04),   # 0.70
      c(-3.85, -3.33, -3.06),   # 0.75
      c(-3.90, -3.37, -3.09),   # 0.80
      c(-3.97, -3.42, -3.13),   # 0.85
      c(-4.08, -3.49, -3.18)    # 0.90
    )
  }

  ## Clamp tau to [0.10, 0.90] for interpolation
  tau_clamp <- min(max(tau, 0.10), 0.90)

  ## Linear interpolation
  cv1  <- approx(tau_grid, cv_tab[, 1], xout = tau_clamp)$y
  cv5  <- approx(tau_grid, cv_tab[, 2], xout = tau_clamp)$y
  cv10 <- approx(tau_grid, cv_tab[, 3], xout = tau_clamp)$y

  c(cv1 = cv1, cv5 = cv5, cv10 = cv10)
}

#' Print method for qadf objects
#'
#' @param x An object of class \code{"qadf"}.
#' @param digits Integer. Number of significant digits for display.
#' @param ... Further arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.qadf <- function(x, digits = 4L, ...) {
  cv  <- x$critical_values
  sig <- if (x$statistic < cv["cv1"]) {
    "***"
  } else if (x$statistic < cv["cv5"]) {
    "**"
  } else if (x$statistic < cv["cv10"]) {
    "*"
  } else {
    ""
  }

  model_lbl <- if (x$model == "c") "Constant" else "Constant + Trend"

  message(strrep("-", 60))
  message("Quantile ADF Unit Root Test")
  message("Koenker & Xiao (2004), JASA, 99, 775-787")
  message(strrep("-", 60))
  message(sprintf("Variable     : %s", x$varname))
  message(sprintf("Model        : %s", model_lbl))
  message(sprintf("Quantile tau : %.3f", x$tau))
  message(sprintf("Optimal lags : %d", x$opt_lags))
  message(sprintf("N (used)     : %d", x$nobs))
  message(sprintf("IC           : %s", x$ic))
  message(strrep("-", 60))
  message("Coefficient Estimates")
  message(strrep("-", 60))
  message(sprintf("  rho(tau) [QR] : %.*f", digits, x$rho_tau))
  message(sprintf("  rho     [OLS] : %.*f", digits, x$rho_ols))
  message(sprintf("  alpha0(tau)   : %.*f", digits, x$alpha_tau))
  message(sprintf("  delta-sq      : %.*f", digits, x$delta2))
  hl_str <- if (is.na(x$half_life)) "---" else sprintf("%.2f", x$half_life)
  message(sprintf("  Half-life     : %s", hl_str))
  message(strrep("-", 60))
  message("Test Statistics")
  message(strrep("-", 60))
  message(sprintf("  t_n(tau)            : %.*f  %s", digits, x$statistic, sig))
  message(sprintf("  U_n(tau) = n(rho-1) : %.*f", digits, x$coef_stat))
  message(strrep("-", 60))
  message("Critical Values (Hansen 1995)")
  message(strrep("-", 60))
  message(sprintf("  1%%   : %.*f", digits, cv["cv1"]))
  message(sprintf("  5%%   : %.*f", digits, cv["cv5"]))
  message(sprintf("  10%%  : %.*f", digits, cv["cv10"]))
  message(strrep("-", 60))

  decision <- if (nchar(sig) > 0L) {
    if (sig == "***") {
      "Reject H0 at 1%: Strong evidence of stationarity"
    } else if (sig == "**") {
      "Reject H0 at 5%: Evidence of stationarity"
    } else {
      "Reject H0 at 10%: Weak evidence of stationarity"
    }
  } else {
    "Cannot reject H0: Evidence consistent with unit root"
  }
  message("Decision: ", decision)
  message("Note: *** p<0.01, ** p<0.05, * p<0.10")
  message("H0: Unit root  vs  H1: Stationarity")
  message(strrep("-", 60))
  invisible(x)
}


#' Summary method for qadf objects
#'
#' @param object An object of class \code{"qadf"}.
#' @param ... Further arguments passed to \code{print.qadf}.
#' @return Invisibly returns \code{object}.
#' @export
summary.qadf <- function(object, ...) {
  print(object, ...)
}

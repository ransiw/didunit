#' A function to get weighted quantiles
#'
#' @param x a numeric vector of draws
#' @param w a numeric vector of weights the same length as draws
#' @param alpha the miscoverage level.
#' @param lower if TRUE, the lower quantile is reported. Default is TRUE.
#'
#' @return a numeric estimate of the lower or upper quantile
#' @export
#' @examples
#' weighted_quantile(rnorm(100,0,1),alpha=0.05)
#'
#'
weighted_quantile <- function(x, w = NULL, alpha, lower = TRUE) {
  if (is.null(w)) w <- rep(1, length(x))
  if (length(w) != length(x)) stop("x and w must be the same length")
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) stop("alpha must be in [0, 1]")
  if (all(is.na(x))) stop("all x missing")

  # keep entries where x is not missing; require corresponding w to be present
  valid <- !is.na(x)
  if (all(is.na(w[valid]))) stop("w must be non-missing wherever x is non-missing")
  x <- x[valid]
  w <- w[valid]

  if (length(x) == 0) return(NA_real_)
  if (any(w < 0)) stop("weights must be non-negative")
  sw <- sum(w)
  if (!(is.finite(sw)) || sw <= 0) stop("sum of weights must be positive and finite")

  # finite-sample correction (jackknife+-style) based on count of non-missing x
  n <- length(x)
  alpha_tilde <- ceiling((n + 1) * alpha) / (n + 1)
  q_level <- if (lower) alpha_tilde else (1 - alpha_tilde)

  # type-1 weighted quantile: smallest x with cumulative weight >= q_level
  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- w[ord]
  cw <- cumsum(w_sorted) / sw

  # handle exact 0 or 1 gracefully
  if (q_level <= 0) return(x_sorted[1])
  if (q_level >= 1) return(x_sorted[length(x_sorted)])

  idx <- which(cw >= q_level)[1]
  x_sorted[idx]
}



#' A function to get weighted standard deviation
#'
#' @param x a numeric vector of draws
#' @param w a numeric vector of weights the same length as draws
#' @param na.rm if TRUE, ignores missing values. Default is TRUE.
#' @param unbiased if TRUE, a Bessel-type correction is applied. Default is TRUE.
#'
#' @return a numeric estimate of standard deviation
#' @export
#' @examples
#' weighted_sd(rnorm(100,0,1))
weighted_sd <- function(x, w = NULL, na.rm = FALSE, unbiased = TRUE) {
  if (is.null(w)) w <- rep(1, length(x))
  if (length(w) != length(x)) stop("length(w) must equal length(x)")
  if (sum(!is.na(x)) == 0) return(NA_real_)

  if (na.rm) {
    # keep entries where x is not missing; require corresponding w to be present
    valid <- !is.na(x)
    if (all(is.na(w[valid]))) stop("w must be non-missing wherever x is non-missing")
    x <- x[valid]
    w <- w[valid]
  } else if (any(is.na(x) | is.na(w))) {
    return(NA_real_)
  }

  if (any(w < 0)) stop("weights must be non-negative")

  W  <- sum(w)
  if (W == 0) return(NA_real_)
  mu <- sum(w * x) / W
  v_num <- sum(w * (x - mu)^2)

  if (unbiased) {
    # unbiased (Bessel-type) correction for general weights
    denom <- W - sum(w^2) / W
    if (denom <= 0) return(NA_real_)
    sqrt(v_num / denom)
  } else {
    sqrt(v_num / W)  # population-style
  }
}

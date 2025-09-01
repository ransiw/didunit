#' @title The average treatment effect calculated using doubly-robust regression.
#'
#' @param index a numeric vector that identifies each unit.
#' @param y1 the outcome vector of the post-treatment period.
#' @param y0 the outcome vector of the pre-treatment period.
#' @param D a binary vector indicating treatment.
#' @param covariates a matrix of covariates, which must include covariates.
#' @param i.weights sample weights, not to be confused with inverse propensity weights.
#' @param dist.family the distribution of the outcome variable. Defaults to gaussian.
#' @param alpha the confidence level for the interval estimates. Defaults to 0.05.
#' @param conformal_split the size of the holdout set when cross-fitting. Defaults to jackknife+ (i.e. 1)
#' @param weighted_conformal whether the conformal scores should be weighted
#'
#' @return a list containing the calculated estimate `ATT`, the conformal interval `lci` and `uci` and a matrix containing the conformal scores and inverse propensity weights.
#' @export
#'
#' @examples
#' # A dataset with a single treated unit
#' data <- sim_data(10,10,1,1, cohort = rep(1,21), pretreat = 1, posttreat = 0)
#' datapre <- data[data$time==9,]
#' data <- data[data$time==10,]
#' # Run the function
#' drobject <- drreg_att(index=data$unit,y1=data$y,y0=datapre$y,D=data$posttreat, covariates = NULL)
#' # Report the ATT
#' drobject$ATT
drreg_att <- function(index,y1, y0, D, covariates, i.weights = NULL, dist.family = stats::gaussian(link = "identity"), alpha = 0.05, conformal_split = 1, weighted_conformal = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(index)
  # generate deltaY
  deltaY <- as.vector(y1 - y0)
  # Covariate vector
  if(is.null(covariates)){
    int.cov <- as.matrix(rep(1,n))
  } else{
    int.cov <- as.matrix(covariates)
  }
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  # Normalize weights
  i.weights <- i.weights/sum(i.weights)

  # Calculate outcome model
  control_filter <- (D == 0)
  reg.coeff <- stats::coef(fastglm::fastglm(
    x = int.cov[control_filter, , drop = FALSE],
    y = deltaY[control_filter],
    weights = i.weights[control_filter],
    family = dist.family
  ))
  out.delta <- as.vector(tcrossprod(reg.coeff, int.cov))

  # Calculate the propensity score model
  PS <- suppressWarnings(fastglm::fastglm(
    x = int.cov,
    y = D,
    family = stats::binomial(),
    weights = i.weights,
    intercept = FALSE,
    method = 3
  ))
  if(PS$converged == FALSE){
    warning("Propensity score estimation did not converge.")
  }
  ps.fit <- stats::fitted(PS)
  # Do not divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-6)
  ipw <- ps.fit / (1 - ps.fit)

  # Calculate the att
  w.treat <- i.weights * D
  w.cont <- i.weights * (1 - D) * ipw

  dr.att.treat <- w.treat * (deltaY - out.delta)
  dr.att.cont <- w.cont * (deltaY - out.delta)

  eta.treat <- mean(dr.att.treat) / mean(w.treat)
  eta.cont <- mean(dr.att.cont) / mean(w.cont)

  dr.att <-   eta.treat - eta.cont

  # Create the vector to save the conformal scores
  CS_full <- rep(NA, length(index))
  M_full <- rep(NA, length(index))
  ipw_full <- rep(NA, length(index))

  # Split the index to holdout
  split_index <- rep(0,length(index))
  split_index[which(D==0)] <- sample(1:floor(length(index[D==0])/conformal_split),length(index[D==0]),replace = TRUE)

  # Populate the matrices
  for (i in unique(split_index)){
    if (i==0){
      next
    } else{
      Ntrain = index[split_index != i & D==0]
      if (length(Ntrain)<dim(int.cov)[2]){
        next
      } else{
        control_filter = index %in% Ntrain

        reg.coeff <- stats::coef(fastglm::fastglm(
          x = int.cov[control_filter, , drop = FALSE],
          y = deltaY[control_filter],
          weights = i.weights[control_filter],
          family = dist.family
        ))
        out.delta <- as.vector(tcrossprod(reg.coeff, int.cov))

        # Calculate the propensity score model
        PS <- suppressWarnings(fastglm::fastglm(
          x = rbind(int.cov[D==1, , drop = FALSE],int.cov[control_filter, , drop = FALSE]),
          y = c(D[D==1],D[control_filter]),
          family = stats::binomial(),
          weights = c(i.weights[D==1],i.weights[control_filter]),
          intercept = FALSE,
          method = 3
        ))
        if(PS$converged == FALSE){
          warning("Propensity score estimation did not converge.")
        }
        ps.fit <- stats::fitted(PS)
        # Do not divide by zero
        ps.fit <- pmin(ps.fit, 1 - 1e-6)
        ipw0 <- ps.fit / (1 - ps.fit)
        ipw0 <- ipw0[-1]

        Ncalib <- index[split_index==i]
        M_full[which(index %in% Ncalib)] <- mean(out.delta[D==1]*i.weights[D==1])/mean(i.weights[D==1])+mean((deltaY-out.delta)[which(index %in% Ntrain)]*i.weights[which(index %in% Ntrain)]*ipw0)/mean(i.weights[which(index %in% Ntrain)]*ipw0)
        CS_full[which(index %in% Ncalib)] <- deltaY[which(index %in% Ncalib)] - out.delta[which(index %in% Ncalib)]
      }
    }
  }

  # Calculate the conformal interval
  if (weighted_conformal){
    ipw_full = ipw
    CSql = weighted_quantile(M_full-abs(CS_full),ipw,alpha = alpha)
    CSqu = weighted_quantile(M_full+abs(CS_full),ipw,alpha = alpha, lower=FALSE)
  } else{
    ipw_full = rep(1,length(index))
    CSql = weighted_quantile(M_full-abs(CS_full),alpha = alpha)
    CSqu = weighted_quantile(M_full+abs(CS_full),alpha = alpha, lower=FALSE)
  }

  CSmat <- data.frame(index = index, CS = CS_full, M = M_full, ipw = ipw_full)

  return(list(ATT = dr.att, lci = mean(deltaY[D==1])-CSqu, uci = mean(deltaY[D==1])-CSql, conf.score.weights = CSmat))

}

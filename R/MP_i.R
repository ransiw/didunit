#' @title MP_i
#'
#' @description Multi-period objects that hold results for unit-time average treatment effects
#'
#' @param id which unit the unit-time treatment effects are for
#' @param group which group (defined by period first treated) an unit-time average treatment effect is for
#' @param t which time period a group-time average treatment effect is for
#' @param att the estimate for unit \code{id} and time
#'  period \code{t}
#' @param lci lower confidence interval for att
#' @param uci upper confidence interval for att
#' @param CS the matrix of conformal scores to produce the aggregate intervals
#' @param IPW the matrix of weights for weighted conformal intervals
#' @param M the matrix of counterfactual predictions
#' @param n the number of unique cross-sectional units (unique values of idname)
#' @param aggite an aggregate treatment effects object
#' @param alp the significance level, default is 0.05
#' @param ipwqual the maximum propensity score
#' @param attcalc similar to att but does not remove estimates when there are propensity score problems
#' @param baseline the baseline of the id
#' @param deltaY the outcome difference from the baseline of the id
#' @param count takes a count of the number of units in the component estimate
#' @param DIDparams a [`DIDparams_i`] object.
#'  of the call to [att_it()].
#'
#' @return MP_i object
#' @export
#'
#' @examples
#' # Helper function for [att_it()]. See documentation of that function for an example.
#'
#'
MP_i <- function(id, group, t, att, lci, uci, CS, IPW, M, n=NULL, aggite=NULL, alp = 0.05, ipwqual=NULL, attcalc=NULL, baseline=NULL, deltaY = NULL, count=NULL, DIDparams=NULL) {
  out <- list(id=id,group=group, t=t, att=att, lci=lci, uci=uci,
              CS = CS, IPW =IPW, M=M, n=n, aggite=aggite, alp = alp, ipwqual=ipwqual,attcalc=attcalc, count=count, baseline=baseline, deltaY = deltaY,
              DIDparams=DIDparams, call=DIDparams$call)
  class(out) <- "MP_i"
  out
}


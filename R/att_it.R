#' Unit-Time Treatment Effects
#'
#' @description
#' Computes the 2x2 DID estimate for each treated unit and each time the treated unit is observed and allows for unit-level covariate matching in a staggered treatment setting.
#'
#'
#' @param yname the name of the outcome variable in the `data`
#' @param tname the name of the column containing the time periods in `data`. Must be numeric
#' @param idname the individual (cross-sectional unit) id name in `data`. Cannot be missing and must be numeric.
#' @param gname the name of the variable in `data` that
#'  contains the first period when a particular observation is treated.
#'  This should be a positive number for all observations in treated groups.
#'  It defines which "treatment-group" a unit belongs to. A zero (0) value is considered to be a unit that is never treated.
#'  A unit treated after the last period in `tname` is also considered never treated.
#' @param cohort the name of blocking variable in `data`, specified if comparisons should only be within pre-treatment cohort.
#'  If blocking is required across the entire trajectory, this must be imposed on the input to `data` prior to running this function.
#' @param xformla a formula for the covariates to include in the
#'  model.  It should be of the form `~ X1 + X2`. An intercept is automatically included. Default
#'  is NULL which is equivalent to `xformla=~1`.
#' @param customnames the names of additional aggregation variables in `data`.
#' @param data the name of the data.frame that contains the data.
#' @param panel if TRUE, the data is coerced to a balanced panel on the. Default is TRUE. No provision for repeated cross-sections differences-in-differences.
#' @param fixedbase if the user needs to fix a base period to be a particular period, use this option. Must be numeric. Default is NULL.
#' @param nobase if TRUE, the baseline vectors are set to zero. Estimates are not differences-in-differences. Default is FALSE.
#' @param control_group which units to use the control group.
#'  The default is "nevertreated" which sets the control group
#'  to be the group of units that never participate in the
#'  treatment.  This group does not change across groups or
#'  time periods.  The other option is to set
#'  `group="notyettreated"`.  In this case, the control group
#'  is set to the group of units that have not yet participated
#'  in the treatment in that time period.  This includes all
#'  never treated units, but it includes additional units that
#'  eventually participate in the treatment, but have not
#'  participated yet.
#' @param anticipation the number of time periods before participating
#'  in the treatment where units can anticipate participating in the
#'  treatment and therefore it can affect their untreated potential outcomes.
#' @param weightsname the name of the column containing the sampling weights.
#'  If not set, all observations have same weight.
#' @param weightfs boolean for whether the att's should be weighted. Default
#' is FALSE so that only the aggite step is weighted unless specified.
#' @param alp the significance level, default is 0.05
#' @param bstrap a deprecated option. Only conformal intervals are reported.
#' @param cband a deprecated option. Does not produce uniform conformal bands.
#' @param biters a deprecated option, as `bstrap`.
#' @param clustervars a deprecated option, as `bstrap`.
#' @param dist.family the assumed distribution for outcome regressions. Accepts any family object of the stats::glm() function, and default is gaussian().
#' @param weighted_conformal option for whether the conformal intervals will be weighted. Default is FALSE.
#' @param conformal_split the size of the holdout set in the cross-fitting procedure. Default is 1 for a leave-one-out jackknife. Increase for faster computation.
#' @param est_method the method to compute group-time average treatment effects.  The default is the Doubly-Robust "dr" Other built-in methods
#' include "ipw" for inverse probability weighting (Hajek type) and "reg" for
#' first step regression estimators.
#' @param overlap the treatment of units that violate overlap conditions when the `est_method` is "dr" or "ipw".
#' The default, "trim", is to drop the treated unit but report the calculated ATT for further analysis. Overlap is violated if the maximum pscore exceeds 0.999
#' The other option, "retain" retains these units for inference.
#' @param base_period whether to use a "varying" base period or a
#'  "universal" base period for placebo tests. A varying base period calculates a pseudo-ATT for every two consecutive pre-treatment periods.
#'  A universal base period fixes the base period to always be (g-anticipation-1).  Either choice results in the same
#'  post-treatment estimates of ATT(i,t)'s.
#' @param print_details whether or not to show details/progress of computations.
#'   Default is `FALSE`.
#' @param pl whether or not to use parallel processing
#' @param cores the number of cores to use for parallel processing
#'
#' @return An [`MP_i`] object containing unit-time treatment effects `att`, the interval estimates `lci` and `uci` and the conformal scores for the aggregation stage.
#' @export
#'
#' @examples
#' attobject <- att_it(yname = "y", tname = "time", gname="treatg", idname="unit", data=sim_data())
#'
att_it <- function(yname,
                   tname,
                   idname,
                   gname,
                   cohort = NULL,
                   xformla = NULL,
                   customnames = NULL,
                   data,
                   panel=TRUE,
                   fixedbase=NULL,
                   nobase=FALSE,
                   control_group=c("nevertreated","notyettreated"),
                   anticipation=0,
                   weightsname=NULL,
                   weightfs=FALSE,
                   alp=0.05,
                   bstrap=TRUE,
                   cband=TRUE,
                   biters=1000,
                   dist.family = stats::gaussian(link = "identity"),
                   weighted_conformal = FALSE,
                   conformal_split = 1,
                   clustervars=NULL,
                   est_method="dr",
                   overlap="trim",
                   base_period="varying",
                   print_details=FALSE,
                   pl=FALSE,
                   cores=1){

  dp <- pre_process_did_i(yname=yname,
                          tname=tname,
                          idname=idname,
                          gname=gname,
                          cohort = cohort,
                          xformla=xformla,
                          customnames = customnames,
                          data=data,
                          panel=panel,
                          fixedbase=fixedbase,
                          nobase=nobase,
                          control_group=control_group,
                          anticipation=anticipation,
                          weightsname=weightsname,
                          weightfs = weightfs,
                          alp=alp,
                          bstrap=bstrap,
                          cband=cband,
                          biters=biters,
                          dist.family = dist.family,
                          weighted_conformal = weighted_conformal,
                          conformal_split = conformal_split,
                          clustervars=clustervars,
                          est_method=est_method,
                          overlap=overlap,
                          base_period=base_period,
                          print_details=print_details,
                          pl=pl,
                          cores=cores,
                          call=match.call()
  )

  results <- compute.att_it(dp)

  # extract ATT(g,t) and influence functions
  attgt.list <- results$attgt.list
  CS <- results$CS
  IPW <- results$IPW
  M <- results$M

  # process results
  attgt.results <- process_attit(attgt.list)
  id <- attgt.results$id
  group <- attgt.results$group
  att <- attgt.results$att
  tt <- attgt.results$tt
  ipwqual <- attgt.results$ipwqual
  lci <- attgt.results$lci
  uci <- attgt.results$uci
  attcalc <- attgt.results$attcalc
  count <- attgt.results$count
  baseline <- attgt.results$baseline
  deltaY <- attgt.results$deltaY

  # number of cross sectional units
  n <- dp$n

  # set the critical value to NULL
  cval <- NULL

  return(MP_i(id=id ,group=group, t=tt, att=att, lci=lci, uci=uci, CS = CS, IPW = IPW, M=M, n=n, alp = alp, ipwqual=ipwqual,attcalc=attcalc,baseline=baseline,deltaY=deltaY,count=count,DIDparams=dp))
}

#' @title Pair Aggregate Unit-Time Average Treatment Effects
#'
#' @description A function to take unit-time average treatment effects
#'  and aggregate them into sub-effects where there is more than one aggregating type.
#'  Possible aggregations include pairs of any of the following:
#'  "group" or a custom aggregation name with "dynamic".
#'
#' @param MP an MP_i object (i.e., the results of the [att_it()] method)
#' @param type Only allows "group" (default) and custom aggregators.
#' @param type2 Only allows "dynamic".
#' @param balance_e If set (and if one computes dynamic effects), it balances
#'  the sample with respect to event time.  For example, if `balance.e=2`,
#'  `aggite` will drop groups that are not exposed to treatment for
#'  at least three periods. (the initial period when `e=0` as well as the
#'  next two periods when `e=1` and the `e=2`).  This ensures that
#'  the composition of groups does not change when event time changes.
#' @param min_e For event studies, this is the smallest event time to compute
#'  dynamic effects for.  By default, `min_e = -Inf` so that effects at
#'  all lengths of exposure are computed.
#' @param max_e For event studies, this is the largest event time to compute
#'  dynamic effects for.  By default, `max_e = Inf` so that effects at
#'  all lengths of exposure are computed.
#' @param na.rm Logical value. If TRUE, removes missing values from analyses. Default is FALSE.
#' @param bstrap Does not apply. To control aggregation assumptions see `indep`.
#' @param indep Logical value. If TRUE, cross-sectional units are aggregated under an independence assumption. Results in a tighter bound. Default is TRUE.
#' @param biters Does not apply.
#' @param cband Does not apply. Set `indep = FALSE` for a conservative bound.
#' @param alp the significance level, default is value set in the MP_i object.
#' @param clustervars Does not apply.
#'
#' @return An [`AGGITEobj`] object that holds the results from the aggregation.
#'
#' @export
#'
#' @examples
#' # first run the att_it() function
#' simdata = sim_data(posttreat=6)
#' attobject = att_it(yname = "y", tname = "time", gname = "treatg", idname ="unit", data = simdata)
#'
#' # aggregate all post-treatment effects of the group and dynamic level
#' agtobject = aggite2(attobject,type="group",type2="dynamic")
#'
#'
#'
aggite2 <- function(MP,
                    type = "group",
                    type2 = "dynamic",
                    balance_e = NULL,
                    min_e = -Inf,
                    max_e = Inf,
                    na.rm = FALSE,
                    indep = TRUE,
                    bstrap = NULL,
                    biters = NULL,
                    cband = NULL,
                    alp = NULL,
                    clustervars = NULL
) {

  call <- match.call()

  compute.aggite2(MP = MP,
                  type = type,
                  type2 = type2,
                  balance_e = balance_e,
                  min_e = min_e,
                  max_e = max_e,
                  na.rm = na.rm,
                  indep = indep,
                  bstrap = bstrap,
                  biters = biters,
                  cband = cband,
                  alp = alp,
                  clustervars = clustervars,
                  call = call)
}

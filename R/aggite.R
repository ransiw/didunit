#' @title Aggregate Unit-Time Average Treatment Effects
#'
#' @description A function to take unit-time average treatment effects
#'  and aggregate them into a smaller number of parameters.  There are
#'  several possible aggregations including "simple", "dynamic", "unit",
#'  "group","calendar, or a custom aggregation name (e.g. "dosage")."
#'
#' @param MP an MP_i object (i.e., the results of the [att_it()] method)
#' @param type Which type of aggregated treatment effect parameter to compute.
#'   One option is "simple" (this just computes a weighted average of all
#'   unit-time average treatment effects).
#'   Other options are "dynamic" (this computes average effects across
#'   different lengths of exposure to the treatment and is similar to an
#'   "event study"; here the overall effect averages the effect of the
#'   treatment across all positive lengths of exposure); "group" (this
#'   is the default option and
#'   computes average treatment effects across different groups; here
#'   the overall effect averages the effect across different groups); and
#'   "calendar" (this computes average treatment effects across different
#'   time periods; here the overall effect averages the effect across each
#'   time period).
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
#' @param na.rm Logical value to remove missing values from MP object's lci or uci. Default is FALSE.
#' @param bstrap Does not apply. To control aggregation assumptions see `indep`.
#' @param indep Logical value. If TRUE, cross-sectional units are aggregated under an independence assumption. Results in a tighter bound. Default is TRUE.
#' @param biters Does not apply.
#' @param cband Does not apply. Set `indep = FALSE` for a conservative bound.
#' @param alp the significance level, default is value set in the MP_i object.
#' @param clustervars Does not apply.
#'
#' @return An [`AGGITEobj`] object that holds the results from the
#'  aggregation
#'
#' @export
#'
#' @examples
#' # first run the att_it() function
#' simdata = sim_data(posttreat=6, pretreat=2)
#' attobject = att_it(yname = "y", tname = "time", gname = "treatg", idname ="unit", data = simdata)
#'
#' # aggregate all post-treatment effects of each unit
#' agtobject = aggite(attobject,type="unit")
#'
#'
aggite <- function(MP,
                   type = "group",
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

  compute.aggite(MP = MP,
                 type = type,
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

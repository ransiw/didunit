#' Aggregate Treatment Effect Parameters Object
#'
#'
#' @inheritParams aggite2
#' @inheritParams compute.aggite2
#' @param overall.att The estimated overall ATT
#' @param overall.se Standard error for overall ATT. Not reported for conformal intervals.
#' @param overall.lci The lower confidence interval of the overall ATT
#' @param overall.uci The upper confidence interval of the overall ATT
#' @param overall.count The number of cross sectional units used to produce the overall ATT
#' @param egt Holds the length of exposure (for dynamic effects), the group, the unit,
#'  cohort, custom aggregator, or the time period (for calendar time effects)
#' @param egt2 a second aggregation type to hold a second-level aggregator, NULL if no such object
#' @param att.egt The ATT specific to egt
#' @param se.egt The standard error specific to egt from the bootstrap. Not reported for conformal intervals.
#' @param lci.egt The lower confidence interval
#' @param uci.egt The upper confidence interval
#' @param count.egt The number of cross-sectional units used to produce the egt.
#' @param crit.val.egt Always NULL because simultaneous option is turned off.
#' @param inf.function The bootstrap draws specific to the egt. Name is carried over from the did package.
#' @param DIDparams A DIDparams_i object
#'
#' @return an AGGITEobj
#' @export
#'
#' @examples
#' # A helper function for [aggite()]. See examples in documentation in that function.
#'
#'
AGGITEobj <- function(overall.att = NULL,
                      overall.se = NULL,
                      overall.lci = NULL,
                      overall.uci = NULL,
                      overall.count = NULL,
                      type = "simple",
                      type2 = NULL,
                      egt = NULL,
                      egt2 = NULL,
                      att.egt = NULL,
                      se.egt = NULL,
                      lci.egt = NULL,
                      uci.egt = NULL,
                      count.egt = NULL,
                      crit.val.egt = NULL,
                      inf.function = NULL,
                      min_e = NULL,
                      max_e = NULL,
                      balance_e = NULL,
                      indep = NULL,
                      call=NULL,
                      DIDparams=NULL) {

  out <- list(overall.att = overall.att,
              overall.se = overall.se,
              overall.lci = overall.lci,
              overall.uci = overall.uci,
              overall.count = overall.count,
              type = type,
              type2 = type2,
              egt = egt,
              egt2 = egt2,
              att.egt = att.egt,
              se.egt = se.egt,
              lci.egt = lci.egt,
              uci.egt = uci.egt,
              count.egt = count.egt,
              crit.val.egt = crit.val.egt,
              inf.function = inf.function,
              min_e = min_e,
              max_e = max_e,
              balance_e = balance_e,
              indep = indep,
              call = call,
              DIDparams = DIDparams)

  class(out)  <- "AGGITEobj"
  out
}


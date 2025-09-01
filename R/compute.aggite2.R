#' Computes pair aggregated treatment effect parameters
#'
#'
#' @inheritParams aggite
#' @param call The function call to aggite2
#'
#' @return [`AGGITEobj`] object
#'
#' @keywords internal
#'
#' @export
#'
#' @examples
#' # This is a helper function for [aggite2()]. See examples in the documentation there.
#'
#'
compute.aggite2 <- function(MP,
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
                            clustervars = NULL,
                            call = NULL) {

  #-----------------------------------------------------------------------------
  # unpack MP_i object
  #-----------------------------------------------------------------------------
  # load parameters
  group <- MP$group
  t <- MP$t
  id <- MP$id
  att <- MP$att
  deltaY <- MP$deltaY
  dp <- MP$DIDparams
  CS <- MP$CS
  M <- MP$M
  IPW <- MP$IPW
  n <- MP$n

  gname <- dp$gname
  data <- as.data.frame(dp$data)
  tname <- dp$tname
  idname <- dp$idname
  cohort <- dp$cohort
  customnames <- dp$customnames


  if(is.null(biters)){
    biters <- dp$biters
  }
  if(is.null(alp)){
    alp <- dp$alp
  }

  tlist <- dp$tlist
  glist <- dp$glist
  idlist <- dp$idlist
  panel <- dp$panel


  if(!(type %in% c("group",customnames))) {
    stop('`type` must be one of c("group", or custom aggregators)')
  }

  if(!(type2 %in% c("dynamic"))) {
    stop('`type2` must "dynamic')
  }

  if(na.rm){
    notna <- (!is.na(att))
    group <- group[notna]
    t <- t[notna]
    id <- id[notna]
    att <- att[notna]
    deltaY <- deltaY[notna]
    CS <- CS[, notna]
    M <- M[,notna]
    IPW <- IPW[, notna]
    #tlist <- sort(unique(t))

    if(type %in% c("unit", "simple", "group", customnames)){
      idlist <- sort(unique(id))
      # Get the units that have some non-missing ATT(g,t) in post-treatmemt periods
      gnotna <- sapply(idlist, function(g) {
        # look at post-treatment periods for group g
        whichg <- which( (id == g) & (group <= t))
        attg <- att[whichg]
        group_select <- !is.na(mean(attg))
        return(group_select)
      })
      gnotna <- idlist[gnotna]
      # indicator for not all post-treatment ATT(g,t) missing
      not_all_na <- id %in% gnotna
      # Re-do the na.rm thing to update the groups
      group <- group[not_all_na]
      t <- t[not_all_na]
      id <- id[not_all_na]
      att <- att[not_all_na]
      deltaY <- deltaY[not_all_na]
      CS <- CS[, not_all_na]
      M <- M[, not_all_na]
      IPW <- IPW[, not_all_na]
      #tlist <- sort(unique(t))
      # redoing the glist here to drop any NA observations
      glist <- unique(data.frame(id,group))$group
      idlist <- sort(unique(id))
    }
  }


  if((na.rm == FALSE) && base::anyNA(att)) stop("Missing values at att_it found. If you want to remove these, set `na.rm = TRUE'.")


  # if the type is a cohort, create cohort variable of the size of ATT(g,t) cohortlist and check that each unit is uniquely mapped to a cohort
  if (type %in% c(customnames)){
    cohortlist <- unique(data[,c(idname,type)])
    idcohort <- data.frame(id = idlist)
    colnames(idcohort) <- idname
    idcohort$save <- rep(1,nrow(idcohort))
    cohortlist <- merge(cohortlist,idcohort, by=idname, sort=FALSE)
    # drop if save is missing
    idcohort <- cohortlist[!is.na(cohortlist$save), , drop = FALSE]
    # find duplicates
    has_multiple_types <- any(duplicated(idcohort$id) | duplicated(idcohort$id, fromLast = TRUE))
    if (has_multiple_types) {
      stop("Some ids belong to multiple cohorts. Consider dropping duplicates")
    }
    # add this to an att_gt
    idcohortatt <- data.frame(id=id)
    colnames(idcohortatt) <- idname
    idcohortatt = merge(idcohortatt,idcohort, by=idname, sort=FALSE)
    ccohort = idcohortatt[,type]
    cohortlist = sort(unique(idcohort[,type]))
  }


  # recover a data-frame with only cross-sectional observations
  if(panel){
    # data from first period
    dta <- data[ data[,tname]==tlist[1], ]
  }else {
    #aggregate data
    dta <- base::suppressWarnings(stats::aggregate(data, list((data[,idname])), mean)[,-1])
  }

  #-----------------------------------------------------------------------------
  # data organization and recoding
  #-----------------------------------------------------------------------------

  # if the na.rm is FALSE the glist for group should be unique
  if (type == "group"){
    glist <- sort(unique(group))
  }

  # do some recoding to make sure time periods are 1 unit apart
  # and then put these back together at the end
  originalt <- t
  originalgroup <- group
  originalglist <- glist
  originaltlist <- tlist
  # In case g's are not part of tlist
  originalgtlist <- sort(unique(c(originaltlist,originalglist)))
  uniquet <- seq(1,length(unique(originalgtlist)))
  # function to switch from "new" t values to  original t values
  t2orig <- function(t) {
    unique(c(originalgtlist,0))[which(c(uniquet,0)==t)]
  }
  # function to switch between "original"
  #  t values and new t values
  orig2t <- function(orig) {
    new_t <- c(uniquet,0)[which(unique(c(originalgtlist,0))==orig)]
    out <- ifelse(length(new_t) == 0, NA, new_t)
    out
  }
  t <- sapply(originalt, orig2t)
  group <- sapply(originalgroup, orig2t)
  glist <- unique(data.frame(id,group))$group
  if (type == "group"){
    glist <- sort(unique(group))
  }
  tlist <- unique(t)
  maxT <- max(t)

  # Set the weights
  weights.ind  <-  dta$.w

  # which group time average treatment effects are post-treatment
  keepers <- which(group <= t & t<= (group + max_e)) ### added second condition to allow for limit on longest period included in att

  # n x 1 vector of group variable
  G <-  unlist(lapply(dta[,gname], orig2t))

  # change the pg to pi to extract the weights for the treated groups
  pg <- dta[dta[,idname] %in% idlist,".w"]

  # length of this is equal to number of treated units or number of units in idlist
  pgg <- pg

  # same but length is equal to the number of ATT(g,t)
  pg <- pg[match(id, idlist)]


  #-----------------------------------------------------------------------------
  # Compute the event-study estimators
  #-----------------------------------------------------------------------------

  if (type2 == "dynamic") {


    # event times
    # this looks at all available event times
    # note: event times can be negative here.
    # note: event time = 0 corresponds to "on impact"
    #eseq <- unique(t-group)
    eseq <- unique(originalt - originalgroup)
    eseq <- eseq[order(eseq)]

    # if the user specifies balance_e, then we are going to
    # drop some event times and some groups; if not, we just
    # keep everything (that is what this variable is for)
    include.balanced.gt <- rep(TRUE, length(originalgroup))

    # if we balance the sample with respect to event time
    if (!is.null(balance_e)) {
      include.balanced.gt <- (t2orig(maxT) - originalgroup >= balance_e)

      eseq <- unique(originalt[include.balanced.gt] - originalgroup[include.balanced.gt])
      eseq <- eseq[order(eseq)]

      eseq <- eseq[ (eseq <= balance_e) & (eseq >= balance_e - t2orig(maxT) + t2orig(1))]

    }

    # only looks at some event times
    eseq <- eseq[ (eseq >= min_e) & (eseq <= max_e) ]

    if (type == "group"){

      # compute atts that are specific to each group and event time
      egtlist = lapply(glist, function(g) {
        lapply(eseq, function(e){list(egt=g,egt2=e)
        })
      })
      egtlist = do.call(c,egtlist)

      dynamic.i.g <- lapply(glist, function(g){
        lapply(eseq,function(e){
          whichi <- which( (group == g) & (originalt - originalgroup == e) & (include.balanced.gt))
          if (length(whichi)==0){
            list(att = NA,lci=NA,uci=NA, count = 0)
          } else{
            whichiN <- id[which( (group == g) & (originalt - originalgroup == e) & (include.balanced.gt))]
            # calculate the aggregate att and deltaYi
            atti <- att[whichi]
            deltaYi <- deltaY[whichi]
            att.i <- sum(atti*(pg[whichi]/sum(pg[whichi])))
            deltaY.i <- sum(deltaYi*(pg[whichi]/sum(pg[whichi])))
            # aggregate predictions and residuals within time blocks
            idx_by_id <- split(seq_along(whichiN),
                               factor(whichiN, levels = unique(whichiN)))
            pgw <- sapply(idx_by_id, function(idx)
              rowMeans(t(as.matrix(pg[whichi]/sum(pg[whichi])))[, idx, drop = FALSE], na.rm = TRUE)
            )
            Mi <- sapply(idx_by_id, function(idx)
              rowMeans(as.matrix(M[,whichi, drop=FALSE])[, idx, drop = FALSE], na.rm = TRUE)
            )
            CSi <- sapply(idx_by_id, function(idx)
              rowMeans(as.matrix(CS[,whichi, drop=FALSE])[, idx, drop = FALSE], na.rm = TRUE)
            )
            IPWi <- sapply(idx_by_id, function(idx)
              rowMeans(as.matrix(IPW[,whichi, drop=FALSE])[, idx, drop = FALSE], na.rm = TRUE)
            )
            if (indep){
              seMi <- sapply(X = seq_len(ncol(Mi[,drop = FALSE])),FUN = function(j) weighted_sd(Mi[,drop = FALSE][, j], IPWi[,drop = FALSE][, j], na.rm = TRUE, unbiased = TRUE) )
              seCSi <- sapply(X = seq_len(ncol(CSi[,drop = FALSE])),FUN = function(j) weighted_sd(CSi[,drop = FALSE][, j], IPWi[,drop = FALSE][, j], na.rm = TRUE, unbiased = TRUE) )
              LC.i <- (deltaY.i-att.i)-stats::qnorm(1-alp/2)*sqrt(sum((seMi^2+seCSi^2)*(pgw/sum(pgw))^2))
              UC.i <- (deltaY.i-att.i)+stats::qnorm(1-alp/2)*sqrt(sum((seMi^2+seCSi^2)*(pgw/sum(pgw))^2))
            } else{
              # find the Bonferroni mean of the prediction intervals
              LCi <- sapply(X = seq_len(ncol(CSi[,drop = FALSE])),FUN = function(j) weighted_quantile(Mi[,drop = FALSE][, j]-abs(CSi[,drop = FALSE])[, j], IPWi[,drop = FALSE][, j],alpha=alp/length(whichiN),lower = TRUE) )
              UCi <- sapply(X = seq_len(ncol(CSi[,drop = FALSE])),FUN = function(j) weighted_quantile(Mi[,drop = FALSE][, j]+abs(CSi[,drop = FALSE])[, j], IPWi[,drop = FALSE][, j],alpha=alp/length(whichiN),lower = FALSE) )
              LC.i <- sum(LCi*(pgw/sum(pgw)))
              UC.i <- sum(UCi*(pgw/sum(pgw)))
            }
            list(att = att.i,lci=deltaY.i-UC.i, uci=deltaY.i-LC.i, count = length(unique(whichiN)))
          }
        })
      })


      dynamic.att.g <- unlist(BMisc::getListElement(do.call(c,dynamic.i.g), "att"))
      dynamic.lci.g <- unlist(BMisc::getListElement(do.call(c,dynamic.i.g), "lci"))
      dynamic.uci.g <- unlist(BMisc::getListElement(do.call(c,dynamic.i.g), "uci"))
      dynamic.count.g <- unlist(BMisc::getListElement(do.call(c,dynamic.i.g), "count"))


      return(AGGITEobj(overall.att=NULL,
                       overall.se=NULL,
                       overall.lci = NULL,
                       overall.uci = NULL,
                       overall.count = NULL,
                       type=type,
                       type2=type2,
                       egt= sapply(unlist(BMisc::getListElement(egtlist, "egt")),t2orig),
                       egt2 = unlist(BMisc::getListElement(egtlist, "egt2")),
                       att.egt=dynamic.att.g,
                       se.egt=NULL,
                       lci.egt=dynamic.lci.g,
                       uci.egt=dynamic.uci.g,
                       count.egt=dynamic.count.g,
                       crit.val.egt=NULL,
                       inf.function = NULL,
                       indep = indep,
                       call=call,
                       min_e=min_e,
                       max_e=max_e,
                       balance_e=balance_e,
                       DIDparams=dp
      ))
    }

    if (type %in% c(customnames)){

      # compute atts that are specific to each group and event time
      egtlist = lapply(cohortlist, function(g) {
        lapply(eseq, function(e){list(egt=g,egt2=e)
        })
      })
      egtlist = do.call(c,egtlist)

      dynamic.i.c <- lapply(cohortlist, function(h){
        lapply(eseq,function(e){
          whichi <- which( (ccohort == h) & (originalt - originalgroup == e) & (include.balanced.gt))
          if (length(whichi)==0){
            list(att = NA,lci=NA,uci=NA, count = 0)
          } else{
            whichiN <- id[which((ccohort == h) & (originalt - originalgroup == e) & (include.balanced.gt))]
            # calculate the aggregate att and deltaYi
            atti <- att[whichi]
            deltaYi <- deltaY[whichi]
            att.i <- sum(atti*(pg[whichi]/sum(pg[whichi])))
            deltaY.i <- sum(deltaYi*(pg[whichi]/sum(pg[whichi])))
            # aggregate predictions and residuals within time blocks
            idx_by_id <- split(seq_along(whichiN),
                               factor(whichiN, levels = unique(whichiN)))
            pgw <- sapply(idx_by_id, function(idx)
              rowMeans(t(as.matrix(pg[whichi]/sum(pg[whichi])))[, idx, drop = FALSE], na.rm = TRUE)
            )
            Mi <- sapply(idx_by_id, function(idx)
              rowMeans(as.matrix(M[,whichi, drop=FALSE])[, idx, drop = FALSE], na.rm = TRUE)
            )
            CSi <- sapply(idx_by_id, function(idx)
              rowMeans(as.matrix(CS[,whichi, drop=FALSE])[, idx, drop = FALSE], na.rm = TRUE)
            )
            IPWi <- sapply(idx_by_id, function(idx)
              rowMeans(as.matrix(IPW[,whichi, drop=FALSE])[, idx, drop = FALSE], na.rm = TRUE)
            )
            if (indep){
              seMi <- sapply(X = seq_len(ncol(Mi[,drop = FALSE])),FUN = function(j) weighted_sd(Mi[,drop = FALSE][, j], IPWi[,drop = FALSE][, j], na.rm = TRUE, unbiased = TRUE) )
              seCSi <- sapply(X = seq_len(ncol(CSi[,drop = FALSE])),FUN = function(j) weighted_sd(CSi[,drop = FALSE][, j], IPWi[,drop = FALSE][, j], na.rm = TRUE, unbiased = TRUE) )
              LC.i <- (deltaY.i-att.i)-stats::qnorm(1-alp/2)*sqrt(sum((seMi^2+seCSi^2)*(pgw/sum(pgw))^2))
              UC.i <- (deltaY.i-att.i)+stats::qnorm(1-alp/2)*sqrt(sum((seMi^2+seCSi^2)*(pgw/sum(pgw))^2))
            } else{
              # find the Bonferroni mean of the prediction intervals
              LCi <- sapply(X = seq_len(ncol(CSi[,drop = FALSE])),FUN = function(j) weighted_quantile(Mi[,drop = FALSE][, j]-abs(CSi[,drop = FALSE])[, j], IPWi[,drop = FALSE][, j],alpha=alp/length(whichiN),lower = TRUE) )
              UCi <- sapply(X = seq_len(ncol(CSi[,drop = FALSE])),FUN = function(j) weighted_quantile(Mi[,drop = FALSE][, j]+abs(CSi[,drop = FALSE])[, j], IPWi[,drop = FALSE][, j],alpha=alp/length(whichiN),lower = FALSE) )
              LC.i <- sum(LCi*(pgw/sum(pgw)))
              UC.i <- sum(UCi*(pgw/sum(pgw)))
            }
            list(att = att.i,lci=deltaY.i-UC.i, uci=deltaY.i-LC.i, count = length(unique(whichiN)))
          }
        })
      })


      dynamic.att.c <- unlist(BMisc::getListElement(do.call(c,dynamic.i.c), "att"))
      dynamic.lci.c <- unlist(BMisc::getListElement(do.call(c,dynamic.i.c), "lci"))
      dynamic.uci.c <- unlist(BMisc::getListElement(do.call(c,dynamic.i.c), "uci"))
      dynamic.count.c <- unlist(BMisc::getListElement(do.call(c,dynamic.i.c), "count"))


      return(AGGITEobj(overall.att=NULL,
                       overall.se=NULL,
                       overall.lci=NULL,
                       overall.uci=NULL,
                       overall.count=NULL,
                       type=type,
                       type2=type2,
                       egt= unlist(BMisc::getListElement(egtlist, "egt")),
                       egt2 = unlist(BMisc::getListElement(egtlist, "egt2")),
                       att.egt=dynamic.att.c,
                       se.egt=NULL,
                       lci.egt=dynamic.lci.c,
                       uci.egt=dynamic.uci.c,
                       count.egt=dynamic.count.c,
                       crit.val.egt=NULL,
                       inf.function = NULL,
                       indep = indep,
                       call=call,
                       min_e=min_e,
                       max_e=max_e,
                       balance_e=balance_e,
                       DIDparams=dp
      ))

    }


  }
}

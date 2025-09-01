#' @title Compute Unit-Time Average Treatment Effects
#'
#' @description
#' `compute.att_it` does the main work for computing
#'  multiperiod unit-time average treatment effects while balancing for the effects at each time period
#'
#'
#' @param dp a DIDparams_i object
#'
#' @return a list object with the calculated att_it() and corresponding conformal scores
#' @export
#'
#' @examples
#' # This is a helper function for [att_it()], but useful for debugging if problems occur.
#' # See documentation of that function for examples.
#'
compute.att_it <- function(dp) {

  #-----------------------------------------------------------------------------
  # unpack DIDparams_i
  #-----------------------------------------------------------------------------
  data <- as.data.frame(dp$data)
  yname <- dp$yname
  tname <- dp$tname
  idname <- dp$idname
  cohort <- dp$cohort
  xformla <- dp$xformla
  weightsname <- dp$weightsname
  weightfs <- dp$weightfs
  est_method <- dp$est_method
  dist.family <- dp$dist.family
  weighted_conformal <- dp$weighted_conformal
  conformal_split <- dp$conformal_split
  overlap <- dp$overlap
  base_period <- dp$base_period
  alp <- dp$alp
  panel <- dp$panel
  fixedbase <- dp$fixedbase
  nobase <- dp$nobase
  print_details <- dp$print_details
  control_group <- dp$control_group
  anticipation <- dp$anticipation
  gname <- dp$gname
  n  <- dp$n
  nT <- dp$nT
  nG <- dp$nG
  tlist <- dp$tlist
  glist <- dp$glist
  idlist <- dp$idlist

  #-----------------------------------------------------------------------------
  # main computations
  #-----------------------------------------------------------------------------

  # will populate with all att(g,t)
  attgt.list <- list()

  # place holder in lists
  counter <- 1

  # number of time periods
  tlist.length <- length(tlist)
  tfac <- 0

  if (base_period == "varying") {
    tlist.length <- tlist.length - 1
    tfac <- 1
  }

  # Matrix to populate with LOO/CV predictions
  Mstore <- Matrix::Matrix(data = 0, nrow = n, ncol = nG*(nT-tfac), sparse = TRUE)

  # Matrix to populate with LOO/CV residuals
  CSstore <- Matrix::Matrix(data = 0, nrow = n, ncol = nG*(nT-tfac), sparse = TRUE)

  # Matrix to populate with propensity weights
  IPWstore <- Matrix::Matrix(data = 0, nrow = n, ncol = nG*(nT-tfac), sparse = TRUE)

  # never treated option
  nevertreated <- (control_group[1] == "nevertreated")

  if(nevertreated) {
    data$.C <- 1*(data[,gname] == 0)
  }

  # rename yname to .y
  data$.y <- data[,yname]


  # determining if the weights are turned on or off
  if (weightfs==FALSE){
    data$.w <- rep(1, nrow(data))
  }


  # Create a dataframe with id and gname
  dataids = data.frame(ids=idlist,
                       gs = glist)
  colnames(dataids) = c(idname,gname)

  # loop over groups
  for (g in 1:nG) {

    # Set up .G once
    data$.G <- 1*(data[,idname] == idlist[g])


    # loop over time periods
    for (t in 1:tlist.length) {

      #-----------------------------------------------------------------------------
      # Set pret

      # varying base period
      pret <- t

      # universal base period
      if (base_period == "universal") {
        # use same base period as for post-treatment periods
        pret <- utils::tail(which( (tlist+anticipation) < glist[g]),1)
      }

      # fixed base period
      if (!is.null(fixedbase)){
        # fix the base period to the value in fixedbase
        pret <- 1
      }

      # use "not yet treated as control"
      # that is, never treated + units that are eventually treated,
      # but not treated by the current period (+ anticipation)
      if(!nevertreated) {
        data$.C <- 1 * ((data[,gname] == 0) |
                          ((data[,gname] > (tlist[max(t,pret)+tfac]+anticipation)) &
                             (data[,gname] != glist[g])))
      }


      # check if in post-treatment period
      if ((glist[g]<=tlist[(t+tfac)])) {

        # update pre-period if in post-treatment period to
        # be  period (g-delta-1)
        pret <- utils::tail(which( (tlist+anticipation) < glist[g]),1)

        # print a warning message if there are no pre-treatment period
        if (length(pret) == 0) {
          warning(paste0("There are no pre-treatment periods for the idgroup first treated at ", idlist[g], "\nUnits from this group are dropped"))

          # if there are not pre-treatment periods, code will
          # jump out of this loop
          break
        }
      }


      #-----------------------------------------------------------------------------
      # if we are in period (g-1), normalize results to be equal to 0
      # and break without computing anything
      if (base_period == "universal") {
        if (tlist[pret] == tlist[(t+tfac)]) {
          attgt.list[[counter]] <- list(att=NA, id=idlist[g],  group=glist[g], year=tlist[(t+tfac)], ipwqual=NA, lci=NA, uci=NA, post=0, baseline=NA, deltaY = NA, attcalc = NA, count=0)
          counter <- counter+1
          next
        }
      }

      # similarly for a fixedbase
      if (!is.null(fixedbase)) {
        if (tlist[pret] == tlist[(t+tfac)]) {
          attgt.list[[counter]] <- list(att=NA, id=idlist[g],  group=glist[g], year=tlist[(t+tfac)], ipwqual=NA, lci=NA, uci=NA, post=0, baseline=NA, deltaY = NA, attcalc = NA, count=0)
          counter <- counter+1
          next
        }
      }

      # print the details of which iteration we are on
      if (print_details) {
        cat(paste("current period:", tlist[(t+tfac)]), "\n")
        cat(paste("current idgroup:", idlist[g]), "\n")
        cat(paste("current group:", glist[g]), "\n")
        cat(paste("set pretreatment period to be", tlist[pret]), "\n")
      }


      # post treatment dummy variable
      post.treat <- 1*(glist[g] <= tlist[t+tfac])

      # total number of units (not just included in G or C)
      disdat <- data[data[,tname] == tlist[t+tfac] | data[,tname] == tlist[pret],]

      # label the pre- and the post-
      disdat$.pre <- as.numeric(disdat[,tname] == tlist[pret])


      n0 <- nrow(disdat)

      # if cohort is not NULL keep to the cohort that matters

      # Save the cohort as a number
      if (!is.null(cohort)){
        cohort0 <- unlist(disdat[disdat$.G==1 & disdat$.pre==1 ,cohort])[1]
      }

      if (!is.null(cohort)){
        precohort_ids = disdat[,idname][disdat[,cohort]==cohort0 & disdat$.pre==1]
        disdat <- disdat[disdat[,idname] %in% precohort_ids,]
        n0 <- nrow(disdat)
      }

      # pick up the indices for units that will be used to compute ATT(g,t)
      # these conditions are (1) you are observed in the right period and
      # (2) you are in the right group (it is possible to be observed in
      # the right period but still not be part of the treated or control
      # group in that period here
      rightids <- disdat[,idname][ disdat$.G==1 | disdat$.C==1]

      # rightids should be observed pre-treatment (imposing same composition on both sides)
      table_rightids <- table(rightids)

      rightids <- rightids[rightids %in% names(table_rightids[table_rightids==2])]

      # pick up the relevant ids
      disidx <- (data[,idname] %in% rightids) & ( (data[,tname] == tlist[t+tfac]) | (data[,tname]==tlist[pret]))

      # pick up the data that will be used to compute ATT(i,t)
      disdat <- data[disidx,]

      # drop missing factors
      disdat <- droplevels(disdat)

      # give short names for data in this iteration
      G <- disdat$.G
      C <- disdat$.C
      Y <- disdat[,yname]
      post <- 1*(disdat[,tname] == tlist[t+tfac])
      # num obs. for computing ATT(i,t)
      n1 <- sum(G+C)
      w <- disdat$.w

      #-----------------------------------------------------------------------------
      # checks to make sure that we have enough observations
      skip_this_att_gt <- FALSE
      if ( sum(G*post) == 0 ) {
        message(paste0("No units for id ", idlist[g], " in time period ", tlist[t+tfac]))
        skip_this_att_gt <- TRUE
      }
      if ( sum(G*(1-post)) == 0) {
        message(paste0("No units for id ", idlist[g], " in time period ", tlist[t]))
        skip_this_att_gt <- TRUE
      }


      if (skip_this_att_gt) {
        attgt.list[[counter]] <- list(att=NA, id=idlist[g],  group=glist[g], year=tlist[(t+tfac)], ipwqual=NA, lci=NA, uci=NA, post=post.treat, baseline=NA, deltaY = NA, attcalc = NA, count=0)
        counter <- counter+1
        next
      }

      skip_this_att_gt <- FALSE
      if (sum(C*post) == 0) {
        message(paste0("No available control units for id ", idlist[g], " in time period ", tlist[t+tfac]))
        skip_this_att_gt <- TRUE
      }
      if (sum(C*(1-post)) == 0) {
        message(paste0("No available control units for group ", idlist[g], " in time period ", tlist[t]))
        skip_this_att_gt <- TRUE
      }

      if (skip_this_att_gt) {
        attgt.list[[counter]] <- list(att=NA, id=idlist[g],  group=glist[g], year=tlist[(t+tfac)], ipwqual=NA, lci=NA, uci=NA, post=post.treat, baseline=NA, deltaY=NA, attcalc = NA, count=0)
        counter <- counter+1
        next
      }

      # Now force to a panel because it balances both sides
      disdat <- BMisc::panel2cs2(disdat, yname, idname, tname, balance_panel=FALSE)
      # drop missing factors
      disdat <- droplevels(disdat)


      # save the indices for the conformal scores
      disidx <- data.frame(id = data[,idname],tt = data[,tname],idx = disidx)
      disidx <- disidx[disidx$tt == tlist[t+tfac],]
      allids <- data.frame(id=sort(unique(data[,idname])))
      allids <- merge(allids,disidx,by="id",all.x=TRUE,sort=TRUE)
      allids$idx[is.na(allids$idx)] <- FALSE
      disidx <- allids$idx

      # sort disdat by idname too
      disdat <- disdat[order(with(disdat, get(idname))), ]

      # give short names for data in this iteration
      G <- disdat$.G
      C <- disdat$.C

      # handle pre-treatment universal base period differently
      # we need to do this because panel2cs2 just puts later period
      # in .y1, but if we are in a pre-treatment period with a universal
      # base period, then the "base period" is actually the later period
      Ypre <- if(tlist[(t+tfac)] > tlist[pret]) disdat$.y0 else disdat$.y1
      if (nobase){
        Ypre <- 0*Ypre
      }

      Ypost <- if(tlist[(t+tfac)] > tlist[pret]) disdat$.y1 else disdat$.y0
      w <- disdat$.w

      # index
      ind <- allids$id[allids$idx]

      # matrix of covariates
      covariates <- stats::model.matrix(xformla, data=disdat)

      # if using custom estimation method, skip this part
      custom_est_method <- class(est_method) == "function"
      maxpscore = NA

      pscore_problems_likely <- FALSE
      reg_problems_likely <- FALSE

      if (!custom_est_method) {

        # checks for pscore based methods
        if (est_method %in% c("dr", "ipw")) {
          preliminary_logit <- stats::glm(G ~ covariates -1, family=stats::binomial(link="logit"))
          preliminary_pscores <- stats::predict(preliminary_logit, type="response")
          if (max(preliminary_pscores) >= 0.999) {
            pscore_problems_likely <- TRUE
            warning(paste0("overlap condition violated for ", glist[g], " in time period ", tlist[t+tfac]))
          }
          maxpscore = max(preliminary_pscores)
        }

        if (pscore_problems_likely & overlap=="trim") {
          attgt.list[[counter]] <- list(att=NA, id=idlist[g], group=glist[g], year=tlist[(t+tfac)],ipwqual=maxpscore, lci=NA, uci=NA, post=post.treat, baseline=Ypre[G==1], deltaY=Ypost[G==1]-Ypre[G==1], attcalc=NA, count=length(ind))
        }

        # check if can run regression using control units
        if (est_method %in% c("dr", "reg")) {
          control_covs <- covariates[G==0,,drop=FALSE]
          #if (determinant(t(control_covs)%*%control_covs, logarithm=FALSE)$modulus < .Machine$double.eps) {
          if ( rcond(t(control_covs)%*%control_covs) < .Machine$double.eps) {
            reg_problems_likely <- TRUE
            warning(paste0("Not enough control units for id ", idlist[g], " in time period ", tlist[t+tfac], " to run specified regression"))
          }
        }

        if (reg_problems_likely) {
          attgt.list[[counter]] <- list(att=NA, id=idlist[g], group=glist[g], year=tlist[(t+tfac)],ipwqual=NA, lci=NA, uci=NA, post=post.treat, baseline=Ypre[G==1], deltaY = Ypost[G==1]-Ypre[G==1], attcalc=NA, count=length(ind))
          counter <- counter+1
          next
        }
      }

      # remove if insufficient for calculating conformal splits
      if (((length(ind[G==0]))-conformal_split)<(ncol(as.matrix(covariates))+3) | floor(length(ind[G==0])/conformal_split)==1){
        attgt.list[[counter]] <- list(att=NA, id=idlist[g], group=glist[g], year=tlist[(t+tfac)],ipwqual=NA, lci=NA, uci=NA, post=post.treat, baseline=Ypre[G==1], deltaY=Ypost[G==1]-Ypre[G==1], attcalc=NA, count=length(ind))
        counter <- counter+1
        warning(paste0("Not enough control units for id ", idlist[g], " in time period ", tlist[t+tfac], "for cross fitting. Consider reducing conformal_split parameter."))
        next
      }


      #-----------------------------------------------------------------------------
      # code for actually computing ATT(i,t)
      #-----------------------------------------------------------------------------

      if (inherits(est_method,"function")) {
        # user-specified function
        attgt <- est_method(ind,Ypost,Ypre,G,
                            covariates = covariates,
                            i.weights = w,
                            alpha = alp,
                            dist.family = dist.family,
                            weighted_conformal = weighted_conformal,
                            conformal_split = conformal_split)
      } else if (est_method == "ipw") {
        # inverse-probability weights
        attgt <- ipw_att(ind,Ypost,Ypre,G,
                         covariates = covariates,
                         i.weights = w,
                         alpha = alp,
                         dist.family = dist.family,
                         weighted_conformal = weighted_conformal,
                         conformal_split = conformal_split)
      } else if (est_method == "reg") {
        # regression
        attgt <- oreg_att(ind,Ypost,Ypre,G,
                          covariates = covariates,
                          i.weights = w,
                          alpha = alp,
                          dist.family = dist.family,
                          weighted_conformal = weighted_conformal,
                          conformal_split = conformal_split)
      } else {
        # doubly robust, this is default
        attgt <- drreg_att(ind,Ypost,Ypre,G,
                           covariates = covariates,
                           i.weights = w,
                           alpha = alp,
                           dist.family = dist.family,
                           weighted_conformal = weighted_conformal,
                           conformal_split = conformal_split)
      }

      if (pscore_problems_likely & overlap=="trim") {
        attgt.list[[counter]] <- list(att=NA, id=idlist[g], group=glist[g], year=tlist[(t+tfac)],ipwqual=maxpscore, lci=NA, uci=NA, post=post.treat, baseline=Ypre[G==1], deltaY=Ypost[G==1]-Ypre[G==1], attcalc=attgt$ATT, count=length(ind))
        counter <- counter+1
        next
      }

      # save to CSstore, Mstore, IPWstore
      CS <- rep(NA, n)
      M <- rep(NA, n)
      IPW <- rep(NA, n)

      # If ATT is NaN, replace it with NA
      if(is.nan(attgt$ATT)){
        attgt$ATT <- NA
        attgt.list[[counter]] <- list(
          att = NA, id=idlist[g],  group = glist[g], year = tlist[(t+tfac)],ipwqual=maxpscore, lci= NA, uci= NA, post = post.treat, baseline=Ypre[G==1], deltaY=Ypost[G==1]-Ypre[G==1], attcalc = NA, count=length(ind))
      } else{
        attgt.list[[counter]] <- list(
          att = attgt$ATT, id=idlist[g],  group = glist[g], year = tlist[(t+tfac)],ipwqual=maxpscore, lci=attgt$lci, uci=attgt$uci, post = post.treat, baseline=Ypre[G==1], deltaY=Ypost[G==1]-Ypre[G==1], attcalc = attgt$ATT, count=length(ind))
        CS[disidx] <- attgt$conf.score.weights[,"CS"]
        CSstore[,counter] <- CS
        M[disidx] <- attgt$conf.score.weights[,"M"]
        Mstore[,counter] <- M
        IPW[disidx] <- attgt$conf.score.weights[,"ipw"]
        IPWstore[,counter] <- IPW
      }

      # update counter
      counter <- counter+1
    } # end looping over t
  } # end looping over g

  return(list(attgt.list=attgt.list, CS = CSstore, M=Mstore, IPW = IPWstore))
}

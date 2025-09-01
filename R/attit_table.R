#' Issues tabulated att_it() results
#'
#' @param MP an `MP_i` object that contains the results from an `att_it()` call
#'
#' @return a dataframe with the id, group, time, att estimates, and their confidence intervals.
#' @export
#'
#' @examples
#' simdata = sim_data()
#' attobject = att_it(yname = "y", tname = "time", gname = "treatg", idname ="unit", data = simdata)
#' attdf = attit_table(attobject)
#'
#'
attit_table <- function(MP){
  #idname = MP$DIDparams$idname
  #gname = MP$DIDparams$gname
  #tname = MP$DIDPparams$tname
  results = data.frame(id=MP$id,group=MP$group,t=MP$t,att=MP$att,lci=MP$lci,uci=MP$uci)
  return(results)
}

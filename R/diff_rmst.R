#' Difference of RMST
#'
#' @export
diff_rmst <- function(est){
   rmst <- diff(est[, "rmst"])
   sd   <-  sqrt( sum(est[,"sd"]^2) )
   wb_sd <- sqrt( sum(est[,"wb_sd"]^2))
   bind_rows(data.frame(est), data.frame(group = 9, rmst = rmst, sd = sd, wb_sd = wb_sd))
}

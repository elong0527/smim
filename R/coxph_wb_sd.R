

#' Wild Bootstrap Estimator for Survival Probability
wild_variance <- function(time, status, u_time, imp_time, s_mi, phi, phi_id = 1:nrow(phi), st_con_survival, n_b, tau, seed, validate = FALSE){

  # u_time    <- sort(unique(time))

  n    <- length(time)
  n_t  <- length(u_time)
  n_mi <- ncol(imp_time)

  st_y <- outer(time, u_time, ">=")

  # For RMST estimation
  loc <- which( u_time < tau)
  dur <- c(u_time[loc], tau) - c(0, u_time[loc])

  ## WB part 1
  v1_imp <- array(NA, c(n_mi, n, n_t))
  for(kk in 1:n_mi){
    if(validate){
      v1_imp[kk,,] <- outer(c(imp_time[status == 0, kk], rep(0, sum(status == 1))), u_time, ">=") * (1 - st_y)
    }else{
      v1_imp[kk,,] <- outer(imp_time[, kk] * (status == 0), u_time, ">=") * (1 - st_y)
    }
  }
  v1_mean <- apply(v1_imp, c(2,3), mean)

  est_wb1<-matrix(0,n_b, n_t)
  rmst_wb1<-rep(0,n_b)
  for( bb in 1:n_b ){
    for( kk in 1:n_mi){
      if(! is.null(seed)){set.seed(seed + bb * 10000 + kk)}
      u_wb        <- rnorm(n)

      thismmu     <- colSums( (v1_imp[kk,,] - v1_mean)*(1- st_y)* u_wb)/ n / n_mi

      est_wb1[bb,]<- est_wb1[bb,]+ thismmu                         # Survival
      rmst_wb1[bb] <- rmst_wb1[bb]+ sum(c(0,thismmu[loc]) * dur) # RMST
    }
  }

  ## WB part 2
  v2_wb     <- phi + st_y + (1 - st_y) * st_con_survival - matrix(s_mi, nrow = n, ncol = n_t, byrow = TRUE)

  # v2_wb_fit <- phi[phi_id, ] + st_y + (1 - st_y) * st_con_survival - matrix(s_mi, nrow = n, ncol = n_t, byrow = TRUE)
  # v2_wb_imp <- phi[- phi_id, ]
  # v2_wb     <- rbind(v2_wb_fit, v2_wb_imp)

  n_v2 <- nrow(v2_wb)
  est_wb2 <- matrix(NA, n_b, n_t)
  rmst_wb2 <- rep(0, n_b)
  for(bb in 1:n_b){
    if(! is.null(seed)){set.seed(seed + bb)}
    u_wb <- rnorm(n_v2)
    est_wb2[bb, ] <- colMeans(v2_wb * u_wb)
    rmst_wb2[bb]  <- sum(c(0, est_wb2[bb, loc]) * dur)
  }

  est_wb <-  est_wb2 + est_wb1
  rmst_wb <- rmst_wb2 + rmst_wb1

  surv_wb_sd <- sqrt(apply(est_wb1,2,var) + colMeans( v2_wb^2)/n)
  rmst_wb_sd <- sd(rmst_wb)
  list(surv_wb_sd = surv_wb_sd,
       rmst_wb_sd = rmst_wb_sd,
       v1_mean = v1_mean,
       v2_wb = v2_wb)
}




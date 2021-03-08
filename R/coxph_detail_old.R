#' Delta adjusted Cox Imputation Model
#'
#' @inheritParams coxph_mi
coxph_delta <- function(fit, pattern, delta){

  # pattern <- as.numeric(R + 1)
  # delta <- c(1,1,2)
  #
  # attach(coxph_martingale(fit))

  fit_res <- coxph_martingale(fit)

  res <- with(fit_res,{

  fit_s$pattern <- as.numeric(pattern)
  fit_s$delta   <- delta[pattern]

  # CAR
  st_survival     <- hazard_to_surv(st_hazard)

  # CNAR
  st_delta_hazard <- ( st_hazard * st_y + fit_s$delta * st_hazard * (1 - st_y) )
  st_delta_survival <- hazard_to_surv(st_delta_hazard)

  s_c <- rowSums( st_survival* st_dy)
  s_c[which(s_c==0)] <- 1
  st_delta_con_survival <- (fit_s$pattern > 1) * (1 - st_y) * (st_survival/s_c)^fit_s$delta # need modification for S1u/S1C part

  # CNAR martingal
  g0 <- colMeans((1 - st_y) * st_delta_con_survival * fit_s$risk)

  g_term <- st_hazard * (1 - st_y) * fit_s$delta * (1 - fit_s$status)

  g1 <- apply(fit$x, 2, function(x){
    g1_term1 <- hazard_to_surv(g_term * x)
    colMeans( (1 - st_y) * st_delta_con_survival * (- log(g1_term1)))
  })

  g2 <- apply(fit$x, 2, function(x){
    t_s1 <- colMeans(fit_s$risk * st_y * x)     # S1
    t_e  <- t_s1 / t_s0                                # E = S1 / S0

    g2_term1 <- hazard_to_surv( t(t(g_term) * t_e) )
    colMeans( (1 - st_y) * st_delta_con_survival * (- log(g2_term1)))
  })

  phi_term1  <- sp_score %*% inv_imat %*% t(g2 - g1)
  phi_term2  <- t(apply(g0 / t_s0 * st_dm * (1 - st_y),1,cumsum))
  phi         <- phi_term1 - phi_term2

  list(fit_s = fit_s,
       st_delta_survival = st_delta_survival,
       st_delta_con_survival = st_delta_con_survival,
       phi = phi)

  })

  c(res, fit_res[-1])

}



#'  Multiple Imputation with Rubin's rule and Wild Bootstrap
#' @param fit a Cox model object, i.e., the result of \code{coxph}.
#' @param pattern A integer vector of pattern indicator for each subject
#'              * 1  = observed event time
#'              * 2  = censoring at random
#'              * 3+ = censoring not at random
#'
#' @param delta A numeric vector of adjusted delta in the same length of pattern.
#'              * delta = 1 is for observed event time and censoring at random.
#'              * delta > 1 is for censoring time with hazard is larger than the hazard in CAR situation.
#'              * delta < 1 is for censoring time with hazard is smaller than the hazard in CAR situation.
#'
#' @param n_mi Number of multiple imputation
#' @param n_wb Number of bootstrap replications
#' @param cut_point Cutpoint of RMST
#' @param seed Randomization seed
#' @param validate Validation model, default is false. Only used for testing purpose.
#'
#' @export
coxph_mi <- function(fit, pattern, delta, n_mi, n_wb, cut_point, seed = NULL, validate = FALSE){

  # pattern <- as.numeric(R + 1)
  # delta <- c(1,1,2)
  if( any(delta[1:2] != 1) ){
    stop("The parameter delta must be 1 for pattern = 1 (observed event time) and pattern = 2 (censoring at random)")
  }

  if( length(delta) != length(unique(pattern))){
    stop("Length of delta does not match unique number of patterns")
  }

  fit0 <- coxph_delta(fit, pattern, delta)
  # attach(fit0)

  res <- with(fit0, {
    n   <- nrow(fit_s)
    n_t <- nrow(fit_t)
    n_imp <- sum(fit_s$status == 0)


    # Multiple Imputation
    imp_upper <- t(st_delta_survival)[which(t(st_dy))] # the upper bound for u's

    imp_data <- matrix(NA, nrow = n, ncol = n_mi)
    for(kk in 1:n_mi){
      # if(! is.null(seed)){set.seed(seed + kk)}

      if(validate){
        u <- imp_upper/kk
      }else{
        u <- runif(n, 0, imp_upper)
      }

      temp <- rowSums( st_delta_survival >= u)
      imp_data[, kk] <- ifelse(fit_s$status == 0, fit_t$time[temp], fit_s$time)
    }

    # Estimation
    s_adj_mi <- var_adj_mi <- matrix(NA, nrow = n_t, ncol = n_mi)
    for(kk in 1:n_mi){
      s_adj_mi[, kk]    <- colMeans( outer( imp_data[, kk], fit_t$time, ">="))
      var_adj_mi[, kk]  <- apply( outer( imp_data[, kk], fit_t$time, ">="), 2, var )/n
    }

    s_adj <- rowMeans(s_adj_mi)                                   # Survival Curve
    loc <- which( fit_t$time < cut_point)
    dur <- c(fit_t$time[loc], cut_point) - c(0, fit_t$time[loc])
    rmst_adj <- sum(c(1, s_adj[loc]) * dur)                       # RMST

    ## RMST
    imp_data_cut <- pmin(imp_data, cut_point)
    rmst_mi <- colMeans(imp_data_cut)
    rmst_var_mi <- apply(imp_data_cut, 2, var) / n

    ## Variance using Rubin's rule
    var_s_adj <- (n.mi+1)/(n.mi)*apply(s_adj_mi,1,var) + rowMeans(var_adj_mi)
    var_rmst   <- (n.mi+1)/(n.mi)*var(rmst_mi) + mean(rmst_var_mi)

    ## RMST results - Rubin's rule
    res_mi <- data.frame(
      type       = "MI with Rubin",
      rmst       = rmst_adj,
      se_rmst    = sqrt(var_rmst),
      lower      = rmst_adj - qnorm(0.975) * sqrt(var_rmst),
      upper      = rmst_adj + qnorm(0.975) * sqrt(var_rmst)
    )
    rownames(res_mi) <- NULL

    # Variance using wild bootstrap (WB)

    ## WB part 1
    v1_imp <- array(NA, c(n_mi, n, n_t))
    v1_imp1 <- array(NA, c(n_mi, n, n_t))
    for(kk in 1:n_mi){

      v1_imp[kk,,] <- outer(c(imp_data[fit_s$status == 0, kk], rep(0, sum(fit_s$status == 1))), fit_t$time, ">=") * (1 - st_y)

      # v1_imp[kk,,] <- outer(imp_data[, kk] * (fit_s$status == 0), fit_t$time, ">=") * (1 - st_y)
    }
    v1_mean <- apply(v1_imp, c(2,3), mean)

    est_wb1<-matrix(0,n_wb, n_t)
    rmst_wb1<-rep(0,n_wb)
    for( bb in 1:n_wb ){
      for( kk in 1:n_mi){
        # if(! is.null(seed)){set.seed(seed + bb * 10000 + kk)}
        u_wb        <- rnorm(n)

        thismmu     <- colSums( (v1_imp[kk,,] - v1_mean)*(1- st_y)* u_wb)/n / n_mi

        est_wb1[bb,]<- est_wb1[bb,]+ thismmu
        rmst_wb1[bb] <- rmst_wb1[bb]+ sum(c(0,thismmu[loc]) * dur)
      }
    }

    ## WB part 2
    v2_wb <- phi + st_y + (1 - st_y) * st_delta_con_survival - matrix(s_adj, nrow = n, ncol = n_t, byrow = TRUE)

    est_wb2 <- matrix(NA, n_wb, n_t)
    rmst_wb2 <- rep(0, n_wb)
    for(bb in 1:n_wb){
      # if(! is.null(seed)){set.seed(seed + bb)}
      u_wb <- rnorm(n)
      est_wb2[bb, ] <- colMeans(v2_wb * u_wb)
      rmst_wb2[bb]  <- sum(c(0, est_wb2[bb, loc]) * dur)
    }

    est_wb <-  est_wb2 + est_wb1
    rmst_wb <- rmst_wb2 + rmst_wb1

    ## WB RMST results
    res_wb <- data.frame(
      type       = "Wild Bootstrap",
      rmst       = rmst_adj,
      se_rmst    = sqrt(var(rmst_wb)),
      lower      = rmst_adj - quantile(rmst_wb,0.975,type=5),
      upper      = rmst_adj - quantile(rmst_wb,0.025,type=5)
    )
    rownames(res_wb) <- NULL

    ## Survival
    fit_t$cnar_surv <- s_adj
    fit_t$cnar_surv_mi_se <- sqrt(var_s_adj)
    fit_t$cnar_surv_wb_se <- sqrt(apply(est_wb1,2,var) + colMeans( v2_wb^2)/n)

    if(validate){
      list(fit_t = fit_t,
           res_wb = res_wb,
           res_mi = res_mi,
           imp_data = imp_data,
           est_wb1 = est_wb1,
           est_wb2 = est_wb2,
           rmst_wb1 = rmst_wb1,
           rmst_wb2 = rmst_wb2,
           s_adj_mi = s_adj_mi,
           var_adj_mi = var_adj_mi
      )
    }else{
      list(fit_t = fit_t,
           res_wb = res_wb,
           res_mi = res_mi)
    }
  })

  c(res, fit0[names(fit0) != "fit_t"])

}

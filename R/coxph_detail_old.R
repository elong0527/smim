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


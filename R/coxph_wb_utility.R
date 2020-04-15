#' Transfer from Hazard Function to Survival Function
#'
#' @noRd
hazard_to_surv <- function(st_hazard){
  temp <- 1 - st_hazard       # Why do we need to use 1 - st_hazard ??
  temp[temp < 0] <- 1
  t(apply(temp,1,cumprod))
}

#' Calculate Phi Function for Wild Bootstrap
#'
#' @noRd
get_phi <- function(st_y, status, sub, delta, x, st_dm, st_con_survival,
                    sp_score, inv_imat, risk, st_hazard){

  if( (ncol(x) == 1) & (length(unique(x[,1])) == 1) ){
    phi <- st_dm * 0
  }else{

    # G0
    g0 <- colMeans((1 - st_y) * st_con_survival * risk)

    g_term <- st_hazard * (1 - st_y) * delta * (1 - status)

    t_s0 <- colMeans((risk * st_y)[sub, ]) # S0
    # G2 - G1
    g21 <- apply(x, 2, function(x){
      t_s1 <- colMeans( (risk * st_y * x)[sub, ])  # S1
      t_e  <- ifelse(t_s0 == 0, 0, t_s1 / t_s0)    # E = S1 / S0

      g1_term1 <- hazard_to_surv(g_term * x)
      g2_term1 <- hazard_to_surv( t(t(g_term) * t_e) )

      g1_term  <- (1 - st_y) * st_con_survival * (- log(g2_term1) + log(g1_term1) )

      colMeans(g1_term[sub, ])
    })

    # Phi
    phi_term1  <- sp_score %*% inv_imat %*% t(g21)
    g0_t_s0    <- ifelse(t_s0 == 0, 0, g0 / t_s0)
    phi_term2  <- t(apply(g0_t_s0 * st_dm * (1 - st_y),1,cumsum))
    phi         <- phi_term1 - phi_term2
  }
  phi
}

#' Calculate Components for Wild Bootstrap
#'
#' @noRd
coxph_wb_utility <- function(fit, id, time, status, x, pattern, delta, fit_imp = NULL){

  if(is.null(fit_imp)){
    coxph_wb_utility_simple(fit, id, time, status, x, pattern, delta)
  }else{
    coxph_wb_utility_imp(fit, id, time, status, x, pattern, delta, fit_imp)
  }

}

#' Calculate Components for Wild Bootstrap for Self Imputation
#'
#' @noRd
coxph_wb_utility_simple <- function(fit, id, time, status, x, pattern, delta, wild_boot = TRUE){

  .x  <- x
  .surv <- Surv(time, status)
  .time_grid <- sort(unique(time))

  fit_res <- coxph_martingale(fit, .id = id, .surv = .surv, .x = .x, .time_grid = .time_grid)

  # Results from coxph_martingale

  sub <- fit_res$sub                 # subset of data used in fit
  risk <- fit_res$fit_s$risk         # risk
  st_y <- fit_res$st_y
  t_s0 <- fit_res$t_s0
  sp_score <- fit_res$sp_score
  inv_imat <- fit_res$inv_imat
  st_dm <- fit_res$st_dm
  st_hazard   <- fit_res$st_hazard

  # sample size
  n <- nrow(st_y)

  # Survival Function
  st_survival <- hazard_to_surv(st_hazard)

  # Delta adjusted Survival function (subject by event time)
  st_delta_hazard <- ( st_hazard * st_y + delta * st_hazard * (1 - st_y) )
  st_delta_survival <- hazard_to_surv(st_delta_hazard)

  # Delta adjusted Conditional Survival function (subject by event time)
  s_c <- st_delta_survival[ cbind(1:n, pmax(rowSums(st_y), 1)) ]
  s_c[which(s_c==0)] <- 1

  st_delta_con_survival <- (status == 0) * (1 - st_y) * (st_survival/s_c)^delta

  if(wild_boot){
    phi <- get_phi(st_y, status, sub, delta, x, st_dm, st_delta_con_survival,
                   sp_score, inv_imat, risk, st_hazard)
    phi[! id %in% fit$id, ] <- 0
  }else{
    phi = NULL
  }

  # Export
  list(
    st_delta_survival = st_delta_survival[fit$id, ],            # With same dim in fit
    st_delta_con_survival = st_delta_con_survival[fit$id, ],    # with same dim in fit
    phi = phi)                                          # with same dim in db
}


#' Calculate Components for Wild Bootstrap for External Imputation
#'
#' @noRd
coxph_wb_utility_imp <- function(fit, id, time, status, x, pattern, delta, fit_imp, wild_boot = TRUE){

  u_pattern <- unique(pattern)
  n_pattern <- length(u_pattern)
  id_pattern <- lapply(u_pattern, function(x) id[which(pattern == x)])

  stopifnot(length(fit_imp) == n_pattern)

  # Subjects used in fit and fit_imp
  .id <- lapply(fit_imp, function(.fit) .fit$id)
  sub.id <- unique(sort(c(fit$id, unlist(.id))))

  .x  <- x
  .surv <- Surv(time, status)
  .time_grid <- sort(unique(time))

  # MAR imputation
  fit_res <- coxph_martingale(fit, .id = id, .surv = .surv, .x = .x, .time_grid = .time_grid)

  fit_pattern <- list()
  for(i in 1:n_pattern){
    if( all(fit$id %in% fit_imp[[i]]$id) & all(fit_imp[[i]]$id %in% fit$id) ){
      fit_pattern[[i]] <- fit_res
    }else{
      fit_pattern[[i]] <- coxph_martingale(fit_imp[[i]], .id = id, .surv = .surv, .x = .x, .time_grid = .time_grid)
    }
  }

  # Results from coxph_martingale

  sub <- fit_res$sub                 # subset of data used in fit
  risk <- fit_res$fit_s$risk         # risk
  st_y <- fit_res$st_y
  t_s0 <- fit_res$t_s0
  sp_score <- fit_res$sp_score
  inv_imat <- fit_res$inv_imat
  st_dm <- fit_res$st_dm
  st_hazard   <- fit_res$st_hazard

  # Hazard function after imputation
  st_hazard_pattern <- lapply(1:n_pattern, function(x) fit_pattern[[x]]$st_hazard[id_pattern[[x]], ])
  st_hazard_pattern <- do.call(rbind, st_hazard_pattern)
  st_hazard_pattern <- st_hazard_pattern[order(unlist(id_pattern)), ]

  # Risk after imputation
  risk_pattern <- lapply(1:n_pattern, function(x) fit_pattern[[x]]$fit_s$risk[id_pattern[[x]]])
  risk_pattern <- do.call(c, risk_pattern)
  risk_pattern <- risk_pattern[order(unlist(id_pattern))]

  # sample size
  n <- nrow(st_y)

  # Survival Function (MAR)
  st_survival <- hazard_to_surv(st_hazard)

  # Delta adjusted (MNAR) Survival function (subject by event time)
  st_delta_hazard <- ( st_hazard * st_y + delta * st_hazard_pattern * (1 - st_y) )
  st_delta_survival <- hazard_to_surv(st_delta_hazard)

  # Delta adjusted Conditional Survival function (subject by event time)
  s_c <- st_delta_survival[ cbind(1:n, pmax(rowSums(st_y), 1)) ]
  s_c[which(s_c==0)] <- 1

  st_delta_con_survival <- (status == 0) * (1 - st_y) * (st_survival/s_c)^delta

  if(wild_boot){
    phi <- list()
    for(i in 1:n_pattern){
      phi[[i]] <- get_phi(st_y, status, sub, delta, x, st_dm, st_delta_con_survival,
                          fit_pattern[[i]]$sp_score, fit_pattern[[i]]$inv_imat,
                          fit_pattern[[i]]$fit_s$risk, fit_pattern[[i]]$st_hazard)
    }

    phi_pattern <- lapply(1:n_pattern, function(x) phi[[x]][id_pattern[[x]], ])
    phi_pattern <- do.call(rbind, phi_pattern)
    phi_pattern <- phi_pattern[order(unlist(id_pattern)), ]

    # Set phi to 0 if the subject is not in fit or fit_imp
    phi_pattern[! id %in% sub.id, ] <- 0
  }else{
    phi_pattern <- NULL
  }

  # Export
  list(
    st_delta_survival = st_delta_survival[fit$id, ],            # With same dim in fit
    st_delta_con_survival = st_delta_con_survival[fit$id, ],    # with same dim in fit
    phi     = phi_pattern)                                      # with same dim in db

}

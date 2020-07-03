
# Delta adjusted RMST within group
rmst_delta <- function(time, status, x, group, pattern, delta, tau, n_mi = 10, n_b = 100, seed = NULL, wild_boot = TRUE, validate = FALSE){

  id <- 1:length(time)

  # Design Matrix
  x  <- as.data.frame(x)
  names(x) <- paste0("x", 1:ncol(x))

  # Time
  u_time <- sort(unique(time))

  # Group
  u_group <- sort(unique(group))
  n_group <- length(u_group)

  # Fit Cox model by group
  fit_group <- list()
  mi_s_group = mi_res_group = fit_wb_group = list()
  for(i in 1:n_group){
    .id <- id[group == u_group[i]]

    db_id <- data.frame(time = time[.id], status = status[.id], x[.id, ])
    # fit <- coxph(Surv(time, status) ~ 1, data = db_id, x = TRUE, y = TRUE)
    fit <- coxph(Surv(time, status) ~ ., data = db_id, x = TRUE, y = TRUE)
    fit$id <- .id


    fit_wb <- coxph_wb_utility_simple(fit, id = id, time = time, status = status, x = x, pattern = pattern, delta = delta, wild_boot = wild_boot)

    mi_t <- mi_time(fit$y[,1], fit$y[,2], u_time,
                    fit_wb$st_delta_survival,
                    n_mi = n_mi, seed = seed, validate = validate)
    mi_s <- mi_survival(fit$y[,1], u_time, mi_t)
    mi_res <- mi_rmst(mi_t, tau = tau)

    if(wild_boot){
      wb_var <- wild_variance(fit$y[,1], fit$y[,2], u_time,
                              mi_t, mi_s[,1], fit_wb$phi, phi_id = fit$id, fit_wb$st_delta_con_survival,
                              n_b = n_b, tau = tau, seed = seed, validate = validate)
    }else{
      wb_var <- NULL
    }

    mi_s <- cbind(group = u_group[i], time = u_time, mi_s, wb_sd = wb_var$surv_wb_sd)
    mi_res <- cbind(group = u_group[i], mi_res, wb_sd = wb_var$rmst_wb_sd)

    fit_group[[i]] <- fit
    mi_s_group[[i]] <- mi_s
    mi_res_group[[i]] <- mi_res
    fit_wb_group[[i]] <- fit_wb
  }

  list(surv = mi_s_group,
       rmst = do.call(rbind, mi_res_group),
       fit_wb_group = fit_wb_group)
}


rmst_control <-  function(time, status, x, group, ref_grp = 0, pattern, delta, tau, n_mi = 10, n_b = 100, seed = NULL, wild_boot = TRUE, validate = FALSE){

  id <- 1:length(time)

  # Design Matrix
  x  <- as.data.frame(x)
  names(x) <- paste0("x", 1:ncol(x))

  # Time
  u_time <- sort(unique(time))

  # Group
  u_group <- sort(unique(group))
  n_group <- length(u_group)

  # Fit Cox model by group
  fit_group <- list()
  for(i in 1:n_group){
    .id <- id[group == u_group[i]]

    db_id <- data.frame(time = time[.id], status = status[.id], x[.id, ])
    fit <- coxph(Surv(time, status) ~ ., data = db_id, x = TRUE, y = TRUE)
    fit$id <- .id

    fit_group[[i]] <- fit
  }

  # Imputation
  mi_s_group = mi_res_group = fit_wb_group = list()
  fit_ref <- fit_group[[which(u_group == ref_grp)]]

  for(i in 1:n_group){

    fit <- fit_group[[i]]

    if(u_group[i] == ref_grp){
      fit_wb <- coxph_wb_utility_simple(fit, id = id, time = time, status = status, x = x, pattern = pattern, delta = delta)
    }else{
      fit_wb <- coxph_wb_utility_imp(fit, id = id, time = time, status = status, x = x, pattern = pattern, delta = delta,
                                     fit_imp = list(fit, fit, fit_ref))
    }

    mi_t <- mi_time(fit$y[,1], fit$y[,2], u_time,
                    fit_wb$st_delta_survival,
                    n_mi = n_mi, seed = seed, validate = validate)
    mi_s <- mi_survival(fit$y[,1], u_time, mi_t)
    mi_res <- mi_rmst(mi_t, tau = tau)

    if(wild_boot){
      wb_var <- wild_variance(fit$y[,1], fit$y[,2], u_time,
                              mi_t, mi_s[,1], fit_wb$phi, phi_id = fit$id, fit_wb$st_delta_con_survival,
                              n_b = n_b, tau = tau, seed = seed, validate = validate)
    }else{
      wb_var <- NULL
    }

    mi_s <- cbind(group = u_group[i], time = u_time, mi_s, wb_sd = wb_var$surv_wb_sd)
    mi_res <- cbind(group = u_group[i], mi_res, wb_sd = wb_var$rmst_wb_sd)

    fit_group[[i]] <- fit
    mi_s_group[[i]] <- mi_s
    mi_res_group[[i]] <- mi_res
    fit_wb_group[[i]] <- fit_wb

  }

  list(surv = mi_s_group,
       rmst = do.call(rbind, mi_res_group),
       fit_wb_group = fit_wb_group)
}


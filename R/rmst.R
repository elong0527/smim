
#' Delta adjusted RMST within group
#'
#' @export
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

    mi_t <- matrix(0, nrow = length(fit$y), ncol = n_mi)
    for(j in 1:n_mi){

    fit_mi <- fit


    coef_post <- as.numeric(rmvnorm(1, mean = fit_mi$coefficients, sigma = fit_mi$var))
    names(coef_post) <- names(fit_mi$coefficients)
    fit_mi$coefficients <- coef_post

    fit_mi$linear.predictors <- NULL
    fit_mi$linear.predictors <- as.numeric((fit_mi$x) %*% matrix(fit_mi$coefficients, nrow = 2)) -
                             as.numeric(fit_mi$means %*%  matrix(fit_mi$coefficients, nrow = 2))

    fit_mi$residuals <- sample(fit$residuals)
    names(fit_mi$residuals) <- names(fit$residuals)



    fit_wb <- coxph_wb_utility_simple(fit_mi, id = id, time = time, status = status, x = x, pattern = pattern, delta = delta, wild_boot = wild_boot)

    mi_t[, j] <- mi_time(fit$y[,1], fit$y[,2], u_time,
                         fit_wb$st_delta_survival,
                         n_mi = 1, pattern = pattern[.id], seed = seed + j *10000, validate = validate)

    }


    mi_s <- mi_survival(fit$y[,1], u_time, mi_t)
    mi_res <- mi_rmst(mi_t, tau = tau)

    if(wild_boot){
      fit_wb0 <- coxph_wb_utility_simple(fit, id = id, time = time, status = status, x = x, pattern = pattern, delta = delta, wild_boot = wild_boot)

      wb_var <- wild_variance(fit$y[,1], fit$y[,2], u_time,
                              mi_t, mi_s[,1], fit_wb0$phi, phi_id = fit$id, fit_wb0$st_delta_con_survival,
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

#' Control-based imputation for RMST
#'
#' @export
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

    .id <- id[group == u_group[i]]
    fit <- fit_group[[i]]

    if(u_group[i] == ref_grp){
      fit_wb <- coxph_wb_utility_simple(fit, id = id, time = time, status = status, x = x, pattern = pattern, delta = delta)
    }else{
      fit_wb <- coxph_wb_utility_imp(fit, id = id, time = time, status = status, x = x, pattern = pattern, delta = delta,
                                     fit_imp = list(fit, fit, fit_ref))
    }

    mi_t <- mi_time(fit$y[,1], fit$y[,2], u_time,
                    fit_wb$st_delta_survival,
                    n_mi = n_mi, pattern = pattern[.id], seed = seed, validate = validate)
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


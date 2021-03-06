# Multiple Imputation

# Multiple Imputation for Event Time
#
mi_time <- function(time, status, u_time, st_survival, n_mi, pattern, seed = NULL, validate = FALSE){

  # u_time    <- sort(unique(time))

  n   <- length(time)
  n_t <- length(u_time)

  st_y <- outer(time, u_time, ">=")
  imp_upper <- st_survival[ cbind(1:n, pmax(rowSums(st_y), 1)) ]

  imp_time <- matrix(NA, nrow = n, ncol = n_mi)
  for(kk in 1:n_mi){
    if(! is.null(seed)){set.seed(seed + kk)}

    if(validate){
      u <- imp_upper/kk
    }else{
      u <- runif(n, 0, imp_upper)
    }

    temp <- rowSums( st_survival >= u)

    # imp_time[, kk] <- ifelse(status == 0, u_time[temp], time)
    imp_time[, kk] <- ifelse(status == 0 & pattern == 2, u_time[temp], time)
    #imp_time[, kk] <- time
  }

  imp_time
}

# Multiple Imputation Estimator for Survival Probability
mi_survival <- function(time, u_time, imp_time){

  # u_time    <- sort(unique(time))

  n    <- length(time)
  n_t  <- length(u_time)
  n_mi <- ncol(imp_time)

  s_mi <- var_mi <- matrix(NA, nrow = n_t, ncol = n_mi)
  for(kk in 1:n_mi){
    s_mi[, kk]    <- colMeans( outer( imp_time[, kk], u_time, ">="))
    var_mi[, kk]  <- apply( outer( imp_time[, kk], u_time, ">="), 2, var )/n
  }

  s <- rowMeans(s_mi) # Survival Curve

  s_var  <- (n_mi+1)/(n_mi)*apply(s_mi,1,var) + rowMeans(var_mi)

  cbind(survival = s, sd = sqrt(s_var) )
}

# Multiple Imputation Estimator for RMST
mi_rmst <- function(imp_time, tau){
  n    <- nrow(imp_time)
  n_mi <- ncol(imp_time)

  imp_time_tau <- pmin(imp_time, tau)
  rmst_mi <- colMeans(imp_time_tau)
  var_mi <- apply(imp_time_tau, 2, var) / n

  rmst_var   <- (n_mi+1)/(n_mi)*var(rmst_mi) + mean(var_mi)
  cbind(rmst = mean(rmst_mi), sd = sqrt(rmst_var))
}





library(survival)

seed <- 123

# Data Simulation
n <- 500
n_mi <- 5
n_b <- 100

truncate_time <- 3.25
tau <- 3

t.time <- rexp(n)   # event time
c.time <- runif(n, min = 0.1, max = 5)  # censoring time
c.time <- pmin(c.time, truncate_time)

time <- pmin(t.time, c.time)
status <- as.numeric(t.time < c.time)
pattern <- ifelse(status == 1, 1, ifelse(time < 0.5, 2, 3))

# Validation Code

fit_km <- survfit(Surv(time, status) ~ 1)


n_t <- length(fit_km$time)
hazard <- c(fit_km$cumhaz[1], diff(fit_km$cumhaz) )
f_surv <- stepfun(fit_km$time, c(1, fit_km$surv), right = FALSE)

s_c <- f_surv(time * (1 - status))  # Survival Probability after censoring


st_hazard <- matrix(hazard, nrow = n, ncol = n_t, byrow = TRUE)
st_y  <- outer(time, fit_km$time, ">=")
st_dn <- outer(time, fit_km$time, "==") * status
st_dm <- (st_dn - st_hazard) * st_y

st_surv <- matrix(fit_km$surv, nrow = n, ncol = n_t, byrow = TRUE)
st_con_surv <- (1 / outer(s_c, fit_km$surv, "/")) * (1 - st_y) * (1 - status)

g0 <- colSums(st_con_surv)
s0 <- fit_km$n.risk
st_g0_s0 <- matrix(g0/s0, nrow = n, ncol = n_t, byrow = TRUE)
phi <-  - t(apply(st_g0_s0 * st_dm, 1, cumsum) ) / 2

# Multiple Imputation
u_time <- fit_km$time
imp_upper <- st_surv[ cbind(1:n, pmax(rowSums(st_y), 1)) ]

imp_time <- matrix(NA, nrow = n, ncol = n_mi)
for(kk in 1:n_mi){
  u <- runif(n, 0, imp_upper)
  temp <- rowSums( st_surv >= u)
  imp_time[, kk] <- ifelse(status == 0, u_time[temp], time)
}

imp_time_tau <- pmin(imp_time, tau)
rmst_mi <- colMeans(imp_time_tau)
var_mi <- apply(imp_time_tau, 2, var) / n

rmst_var   <- (n_mi+1)/(n_mi)*var(rmst_mi) + mean(var_mi)


s_mi <- matrix(NA, nrow = n_t, ncol = n_mi)
for(kk in 1:n_mi){
  s_mi[, kk]    <- colMeans( outer( imp_time[, kk], u_time, ">="))
}
s_mi <- rowMeans(s_mi) # Survival Curve

# For RMST estimation
loc <- which( u_time < tau)
dur <- c(u_time[loc], tau) - c(0, u_time[loc])

## WB part 1
v1_mean <- matrix(0, n, n_t)
for(kk in 1:n_mi){
  v1_mean <- v1_mean + outer(imp_time[, kk] * (status == 0), u_time, ">=") * (1 - st_y)
}
v1_mean <- v1_mean / n_mi

rmst_wb1<-rep(0,n_b)
for( bb in 1:n_b ){
  for( kk in 1:n_mi){
    u_wb        <- rnorm(n)

    v1_imp <- outer(imp_time[, kk] * (status == 0), u_time, ">=") * (1 - st_y)
    thismmu     <- colSums( (v1_imp - v1_mean)*(1- st_y)* u_wb)/ n / n_mi

    rmst_wb1[bb] <- rmst_wb1[bb]+ sum(c(0,thismmu[loc]) * dur) # RMST
  }
}

v2_wb <- phi + st_y + (1 - st_y) *  st_con_surv - matrix(s_mi, nrow = n, ncol = n_t, byrow = TRUE)

est_wb2 <- matrix(NA, n_b, n_t)
rmst_wb2 <- rep(0, n_b)
for(bb in 1:n_b){
  if(! is.null(seed)){set.seed(seed + bb)}
  u_wb <- rnorm(n)
  est_wb2[bb, ] <- colSums(v2_wb * u_wb) / n
  rmst_wb2[bb]  <- sum(c(0, est_wb2[bb, loc]) * dur)
}

rmst_wb <- rmst_wb2 + rmst_wb1


res <- cbind(rmst = mean(rmst_mi), sd = sqrt(rmst_var), rmst_wb_sd = sd(rmst_wb))

fit <- coxph(Surv(time, status) ~ rep(1,n), x = TRUE, y = TRUE)
res1 <- rmst_delta(time, status, x = fit$x, group = 1, pattern = rep(1, n), delta = rep(1, n), tau = tau, n_mi = n_mi, n_b = n_b, seed = seed)

# phi <-  - t(apply((1 - status) * (1 - st_y) * st_g0_s0 * st_dm, 1, cumsum) )

# Validation Start
fit <- coxph(Surv(time, status) ~ rep(1,n), x = TRUE, y = TRUE)

martingale_val <- coxph_martingale(fit)

expect_equivalent(st_y, martingale_val$st_y)
expect_equivalent(st_hazard, martingale_val$st_hazard)
expect_equivalent(st_dm, martingale_val$st_dm)

fit$id <- 1:n
wb_val <- coxph_wb_utility_simple(fit = fit,
                                  id = fit$id,
                                  time = as.numeric(time),
                                  status = as.numeric(status),
                                  x = fit$x,
                                  pattern = pattern,
                                  delta = c(1,1,1)[pattern])

expect_equivalent(wb_val$st_delta_survival, st_surv)
expect_equivalent(wb_val$st_delta_con_survival, st_con_surv)
expect_equivalent(wb_val$phi, phi)

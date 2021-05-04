##########################################
#  Simulation Start
##########################################

task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
set.seed(task_id)
seed <- task_id

library(survival)
library(dplyr)
library(smim)
library(mvtnorm)

simu_one_group <- function(n, x, beta, lambda, betaC, lambdaC, tau, arm, LL,
                           beta_imp = beta, lambda_imp = lambda, delta = 1, vt = NULL, vc = NULL){

  # Simulate time


  if(is.null(dim(x))){ x <- matrix( x, ncol = 1) }

  if(is.null((vt))){
    v <- runif(n)
  }else{
    v <- vt
  }

  time  <- (- log(1-v) / (lambda * exp(x %*% beta)))

  # Simulate Censoring time
  if(is.null((vc))){
    v <- runif(n)
  }else{
    v <- vc
  }
  censor <- (- log(1-v) / (lambdaC * exp(x %*% betaC)))
  # censor <- (- log(1-v) / (lambdaC) )
  censor <- pmin(censor, LL)

  status <- ifelse(time < censor, 1, 0)
  pattern <- ifelse(time < censor, 0,
                    ifelse(censor == LL, 1, 2))
  # Observed time
  obs <- pmin(time, censor)

  v <- runif(n)
  time2 <- (- log(1-v) / (delta * lambda_imp * exp(x %*% beta_imp)))
  time <- ifelse(pattern == 2, censor + time2, time)

  data.frame(time, censor, status, obs, arm, x, pattern)

}

n_b  <- 100
tau <- 24
LL  <- 40

delta <- 2

beta0 <- c(-0.55, 0.65)
lambda0 <- 0.03

beta1 <- c(0.24, 0.04)
lambda1 <- 0.03

betaC     <- c(0.25, 0.20)
lambdaC   <- 0.01

# parameters

n <- 500
n_mi <- 20

db0 <- simu_one_group(n = n, x = cbind(rnorm(n), rbinom(n, size = 1, prob = 0.15) ), beta = beta0, lambda = lambda0, betaC = betaC, lambdaC = lambdaC, tau = tau, arm = 0, LL = LL, delta = 1)
db1 <- simu_one_group(n = n, x = cbind(rnorm(n), rbinom(n, size = 1, prob = 0.15) ), beta = beta1, lambda = lambda1, betaC = betaC, lambdaC = lambdaC, tau = tau, arm = 1, LL = LL, delta = delta)

db <- rbind(db0, db1)

# True value
n_true <- 10000
dt0 <- simu_one_group(n = n_true, x = cbind(rnorm(n_true), rbinom(n_true, size = 1, prob = 0.15) ), beta = beta0, lambda = lambda0, betaC = betaC, lambdaC = lambdaC, tau = tau, arm = 0, LL = LL, delta = 1)
dt1 <- simu_one_group(n = n_true, x = cbind(rnorm(n_true), rbinom(n_true, size = 1, prob = 0.15) ), beta = beta1, lambda = lambda1, betaC = betaC, lambdaC = lambdaC, tau = tau, arm = 1, LL = LL, delta = delta)

dt <- rbind(dt0, dt1)
pct <- rbind(table(dt0$pattern), table(dt1$pattern), table(dt$pattern)/2) / n_true * 100

if(i == 1){
  print(pct)
}

dt_rmst <- c(
  mean(pmin(subset(dt, arm == 0)$time, tau)),
  mean(pmin(subset(dt, arm == 1)$time, tau))
)
dt_rmst[3] <- diff(dt_rmst)

# Estimated from rmst_delta
delta1 <- seq(1, 5, by = 1)
# delta1 <- 2
diff <- list()
for(iter in 1:length(delta1)){
  db$group <- db$arm

  db$delta_est <- ifelse(db$pattern == 2 & db$group == 1, delta1[iter], 1)
  # Delta-adjusted estimator
  fit <- rmst_delta(db$obs,
                    db$status,
                    x = db[, c("X1", "X2")],
                    group = db$group,
                    pattern = db$pattern,
                    delta = db$delta_est,
                    tau = tau,
                    n_mi = par0$n_mi,
                    n_b = n_b,
                    seed = seed, wild_boot = TRUE)
  # fit$rmst <- cbind(fit$rmst, wb_sd = 0)

  diff[[iter]] <-  data.frame(diff_rmst(fit$rmst), delta = delta1[iter])
}
diff <- bind_rows(diff)

# RMST from survRM2
fit <- survRM2::rmst2(time = db$obs, status = db$status, arm = db$arm, tau = tau)
est <- bind_rows(fit$RMST.arm0$rmst, fit$RMST.arm1$rmst, fit$unadjusted.result[1, ])
sd  <- c(sqrt(fit$RMST.arm0$rmst.var), sqrt(fit$RMST.arm1$rmst.var),
         sqrt(fit$RMST.arm0$rmst.var + fit$RMST.arm1$rmst.var))
est$sd <- sd
rmst2 <- est %>% mutate(group = c(0,1,9), rmst = Est., sd = sd) %>%
  select(group, rmst, sd)
rmst2$delta <- 88

#Output
res <- data.frame(par0, bind_rows(diff, rmst2))
res_true <- data.frame(par0, group = c(0,1,9), true = dt_rmst, pct = pct)


# Save Simulation Results
filename <- paste0(task_id,".Rdata")
save(res, res_true, file = filename)


#----------------------
# HPC code Submission
#----------------------

# module add R/3.6.3
# cd /SFS/scratch/zhanyilo/smim6
# rm *
# qsub -t 1:1000 ~/runr.sh ~/smim/simulation/simu_actg.R

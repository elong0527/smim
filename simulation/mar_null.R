#---------------------
# Simulation seed
#---------------------

task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))

# Set up Simulation Enviroment

# task_id <- 1
seed <- task_id
###############################

library(purrr)
library(dplyr)
library(miboot)
library(survRM2)
source("simu_data.R")



simu_mar_null <- function(n_true, n_grp, delta_mnar, trunc_time, tau,
                     lambda, beta, lambda_cen, beta_cen, lambda_control = NULL, beta_control = NULL,
                     n_mi, n_b,
                     seed = NULL){

  # n_true = 10000
  # n_grp = c(250)
  # delta_mnar = c(0.8)
  # trunc_time = c(3.25)
  # tau = c(3)
  # lambda = 0.35
  # beta =0.75
  # lambda_cen = 0.15
  # beta_cen = c(0.75)
  # n_mi = 5
  # n_b = 100
  # seed =  1

  # True value
  grp1_true <- simu_mnar_1grp(n_true, delta_mnar, trunc_time, lambda, beta, lambda_cen, beta_cen, lambda_control = NULL, beta_control = NULL, seed = seed)
  grp1_true$group <- 1
  grp0_true <- simu_mnar_1grp(n_true, delta_mnar, trunc_time, lambda, beta, lambda_cen, beta_cen, lambda_control = NULL, beta_control = NULL, seed = seed + 1e4)
  grp0_true$group <- 0

  db_true <- rbind(grp1_true, grp0_true)

  true_value <- db_true %>% group_by(group) %>%
    summarise(rmst_true = mean(pmin(t.time, tau)),
              cen_prop = mean(status == 0),
              pattern_2 = mean(pattern == 2),
              pattern_3 = mean(pattern == 3))

  # Estimation
  grp1 <- simu_mnar_1grp(n_grp, delta_mnar, trunc_time, lambda, beta, lambda_cen, beta_cen, lambda_control = NULL, beta_control = NULL, seed = seed)
  grp1$group <- 1
  grp0 <- grp0_true <- simu_mnar_1grp(n_grp, delta_mnar, trunc_time, lambda, beta, lambda_cen, beta_cen, lambda_control = NULL, beta_control = NULL, seed = seed + 1e4)
  grp0$group <- 0

  db <- rbind(grp1, grp0)

  ## RMST within group
  tmp <- with(db, rmst_delta(time, status, x = rep(1, length(time)), group, pattern, delta, tau, n_mi = n_mi, n_b = n_b, seed = seed))
  tmp$rmst
  # tmp <- list()
  # tmp$rmst <- structure(c(0, 1, 0.947450795954468, 0.895520988424821, 0.0386727153576288,
  #                         0.0352587959938403, 0.0290456519730724, 0.0490021896154576),
  #                       .Dim = c(2L, 4L), .Dimnames = list(NULL, c("group", "rmst", "sd", "wb_sd")))

  ## RMST between group
  diff_rmst <- function(rmst, sd){
    diff <- diff(rmst)
    diff_sd <- sqrt( sum(sd^2) )
    p_val <- 2* (1 - pnorm( abs(diff/diff_sd) ))
    c(diff = diff, sd = diff_sd, p = p_val)
  }

  ## Prepare Simulation Results
  res <- tmp$rmst

  diff_res <- rbind( diff_rmst(tmp$rmst[,"rmst"], tmp$rmst[, "sd"]),
                     diff_rmst(tmp$rmst[,"rmst"], tmp$rmst[, "wb_sd"]))
  diff_res <- data.frame(type = c("rubin", "wb"), diff_res)

  list(res = res, diff_res = diff_res, true_value = true_value)

}

par <- expand.grid(
  n_true = 10000,
  n_grp = c(250),
  delta_mnar = c(0.1, 1, 2),
  trunc_time = c(5.25),
  tau = c(5),
  lambda = 0.35,
  beta =0.75,
  lambda_cen = 0.35,
  beta_cen = c(0.75),
  n_mi = c(10),
  n_b = 100,
  seed =  task_id
)

# par <- expand.grid(
#   n_true = 10000,
#   n_grp = c(250 , 500),
#   delta_mnar = c(0.8, 1, 1.2, 1.5),
#   trunc_time = c(3.25),
#   tau = c(1, 2, 3),
#   lambda = 0.35,
#   beta =0.75,
#   lambda_cen = 0.15,
#   beta_cen = c(0.75, 1),
#   n_mi = c(5, 10),
#   n_b = 100,
#   seed =  task_id
# )

par$result <- pmap(par, simu_mar_null)

save(par, file = paste0(task_id, ".Rdata") )

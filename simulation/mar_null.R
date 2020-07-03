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
library(smim)
library(survRM2)
source("simu_data.R")



simu_mar_null <- function(n_true, n_grp, delta_mnar, trunc_time, tau,
                     lambda, beta, lambda_cen, beta_cen, lambda_control = NULL, beta_control = NULL,
                     n_mi, n_b,
                     seed = NULL){

  # n_true = 10000
  # n_grp = c(100)
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
  tmp <- with(db, rmst_delta(time, status, x = x, group, pattern, delta, tau, n_mi = n_mi, n_b = n_b, seed = seed))
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



# n <- 1000
# n_mi <- 20   # number of MI
# n_b <- 100   # number of bootstrap
#
# t.time <- rexp(n)   # event time
# # t.time <- rweibull(n, shape = 0.8, scale = 1.5)
# c.time <- runif(n, min = 0.1, max = 5)  # censoring time
#
# time <- pmin(t.time, c.time)
#
# status <- as.numeric(t.time < c.time)
#
# pattern <- ifelse(status == 1, 1, ifelse(time < 0.5, 2, 3))
# group <- rep(c(0,1), each = 500)
#
# #plot(survfit(Surv(time, status) ~ group), col = 1:2)
# #table(status, pattern, time < 0.5)
#
# # MAR analysis
# tau <- 3
#
# ## survRM2 results
# rm2 <- rmst2(time, status, group, tau = tau)
# rm2_res <- data.frame(group = c(0,1), rbind(rm2$RMST.arm0$rmst, rm2$RMST.arm1$rmst))
# rm2_diff <- rm2$unadjusted.result[1,]
#
# ## RMST within group
# delta   <- c(1,1,1)[pattern]  # the third number control delta adjustment for MNAR
# tmp <- rmst_delta(time, status, x = rep(1, length(time)), group, pattern, delta, tau, n_mi = n_mi, n_b = n_b, seed = seed)
# tmp$rmst
#
#
# ## RMST between group
# diff_rmst <- function(rmst, sd){
#   diff <- diff(rmst)
#   diff_sd <- sqrt( sum(sd^2) )
#   p_val <- 2* (1 - pnorm( abs(diff/diff_sd) ))
#   c(diff = diff, sd = diff_sd, p = p_val)
# }
#
# ## Prepare Simulation Results
# res <- tmp$rmst
#
# diff_res <- rbind( diff_rmst(tmp$rmst[,"rmst"], tmp$rmst[, "sd"]),
#                    diff_rmst(tmp$rmst[,"rmst"], tmp$rmst[, "wb_sd"]),
#                    c(rm2_diff["Est."], diff(rm2_diff[2:3])/2/qnorm(0.975), rm2_diff[4] ) )
# diff_res <- data.frame(type = c("rubin", "wb", "survRM2"), diff_res)
#
# save(res, rm2_res, diff_res, file = paste0(task_id, ".Rdata") )


library(purrr)
library(dplyr)
library(tidyr)

path <- "/SFS/scratch/zhanyilo/mi_mar_null/"


res = list()
for(i in 1:1000){
  try({
    load(file.path(path, paste0(i, ".Rdata")))
    res[[i]] <- par
  })
}

library(purrr)
library(dplyr)
library(tidyr)

path <- "/SFS/scratch/zhanyilo/mi_mar_null/"


res = list()
for(i in 1:1000){
  try({
    load(file.path(path, paste0(i, ".Rdata")))
    res[[i]] <- par
  })
}

# db <- do.call(rbind, res)
#
# ## True Value
# res_true <- db %>% mutate(
#   res = map(.$result, "true_value")) %>%
#   select(- result) %>% unnest(res) %>%
#   group_by(n_true, n_grp, delta_mnar, trunc_time, tau, lambda, beta, lambda_cen, beta_cen, n_mi, n_b) %>%
#   summarise(rmst_true = mean(rmst_true))
#
# ## Within Group
# res_group <- db %>% mutate(
#   res = map(.$result, function(x) as.data.frame(x$res))) %>%
#   select(- result) %>% unnest(res) %>%
#   left_join(res_true) %>%
#   group_by(n_true, n_grp, delta_mnar, trunc_time, tau, lambda, beta, lambda_cen, beta_cen, n_mi, n_b, group) %>%
#   summarise(n = n(),
#             rmst_true = mean(rmst_true),
#             rmst_est  = mean(rmst),
#             sd_empirical = sd(rmst),
#             sd_rubin = mean(sd),
#             sd_wb = mean(wb_sd),
#             covr_rubin = mean( rmst_true > (rmst - 1.96 * sd) & rmst_true < (rmst + 1.96 * sd) ),
#             covr_wb =  mean( rmst_true > (rmst - 1.96 * wb_sd) & rmst_true < (rmst + 1.96 * wb_sd) ) )
#
# ## Between Group
# res_diff <- db %>% mutate(
#   res = map(.$result, "diff_res")) %>%
#   select(- result) %>% unnest(res) %>%
#   group_by(n_true, n_grp, delta_mnar, trunc_time, tau, lambda, beta, lambda_cen, beta_cen, n_mi, n_b, type) %>%
#   summarise(n = n(),
#             diff_est = mean(diff),
#             sd_empirical = sd(diff),
#             sd_est = mean(sd),
#             type1 = mean(p < 0.05))
#
#
#
# library(ggplot2)
# qplot(tau, type1, group = type, color = type, data = res_diff) + facet_grid(n_grp ~ n_mi)
#
# res_group_long <- res_group %>% pivot_longer(covr_rubin:covr_wb)
#
# qplot(tau, value, group = name, color = paste(name, group), data = res_group_long) + facet_grid(n_grp + dist ~ n_mi)
#
# save.image(file = "simulation/mcar.Rdata")

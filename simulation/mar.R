#---------------------
# Simulation seed
#---------------------

i.par <- commandArgs(T)
i.par <- as.numeric(i.par)
root.file <- paste0(i.par)

seed <- i.par
###############################


library(miboot)
library(survRM2)

# * Pattern
#   + 1: observed event time
#   + 2: MAR imputation
#   + 3: MNAR imputation
#
# * Delta
#   + 1: no adjustment
#   + larger than 1: hazard is larger than MAR
#   + less than 1: hazard is smaller than MAR

n <- 1000
t.time <- rexp(n)   # event time
c.time <- runif(n, min = 0.1, max = 5)  # censoring time

time <- pmin(t.time, c.time)

status <- as.numeric(t.time < c.time)

pattern <- ifelse(status == 1, 1, ifelse(time < 0.5, 2, 3))
group <- rep(c(0,1), each = 500)

#plot(survfit(Surv(time, status) ~ group), col = 1:2)
#table(status, pattern, time < 0.5)

# MAR analysis
tau <- 3

## survRM2 results
rm2 <- rmst2(time*10, status, group, tau = tau)
rm2_res <- data.frame(group = c(0,1), rbind(rm2$RMST.arm0$rmst, rm2$RMST.arm1$rmst))
rm2_diff <- rm2$unadjusted.result[1,]

## RMST within group
delta   <- c(1,1,1)[pattern]  # the third number control delta adjustment for MNAR
tmp <- rmst_delta(time, status, x = rep(1, length(time)), group, pattern, delta, tau, n_mi = 10, n_b = 100, seed = seed)
tmp$rmst


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
                   diff_rmst(tmp$rmst[,"rmst"], tmp$rmst[, "wb_sd"]),
                   c(rm2_diff["Est."], diff(rm2_diff[2:3])/2/qnorm(0.975), rm2_diff[4] ) )
diff_res <- data.frame(type = c("rubin", "wb", "survRM2"), diff_res)

save(res, rm2_res, diff_res, file = file.path(i.par, ".Rdata") )





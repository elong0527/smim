n <- 1000
t.time <- rexp(n)   # event time
c.time <- runif(n, min = 0.1, max = 5)  # censoring time

time <- pmin(t.time, c.time)

status <- as.numeric(t.time < c.time)

pattern <- ifelse(status == 1, 1, ifelse(time < 0.5, 2, 3))
delta <- c(1,1,1)[pattern]
x <- matrix(1, nrow = n, ncol = 1)

group <- rep(c(0,1), each = 500)

id <- 1:n
u_group <- sort(unique(group))
u_time <- sort(unique(time))

for(i in 1:length(u_group)){
  .id <- id[group == u_group[i]]


  db_id <- data.frame(time = time[.id], status = status[.id], x[.id, ])
  # fit <- coxph(Surv(time, status) ~ 1, data = db_id, x = TRUE, y = TRUE)
  fit <- coxph(Surv(time, status) ~ ., data = db_id, x = TRUE, y = TRUE)
  fit$id <- .id


  fit_wb <- coxph_wb_utility_simple(fit, id = id, time = time, status = status, x = x, pattern = pattern, delta = delta, wild_boot = TRUE)
  fit_wb1 <- coxph_wb_utility_simple(fit, id = fit$id, time = fit$y[,1], status = fit$y[,2], x = fit$x, pattern = pattern[.id], delta = delta[.id], wild_boot = TRUE)

  test_that(paste("Validate expand subject and time in Group", u_group[i]), {
    expect_equivalent(fit_wb1$st_delta_survival, fit_wb$st_delta_survival[,u_time %in% sort(unique(time[.id]))])
    expect_equivalent(fit_wb1$st_delta_con_survival, fit_wb$st_delta_con_survival[,u_time %in% sort(unique(time[.id]))])
    expect_equivalent(fit_wb1$phi, fit_wb$phi[id %in% .id, u_time %in% sort(unique(time[.id])) ] )
  })

}

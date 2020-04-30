
# simu_mcar_1grp(100, 1, 4, 5)

simu_mcar_1grp <- function(n, mean, trunc_time, cen_max, dist = "exp", seed = NULL){

  if(! is.null(seed)){
    set.seed(seed)
  }

  if(dist == "exp"){
    t.time <- rexp(n)   # event time
  }

  if(dist == "weibull"){
    shape <- 0.8
    scale <- mean / gamma(shape)
    t.time <- rweibull(n, shape = shape, scale = scale)
  }

  c.time <- runif(n, min = 0.1, max = cen_max)  # censoring time
  c.time <- pmin(trunc_time, c.time)

  time <- pmin(t.time, c.time)

  status <- as.numeric(t.time < c.time)

  pattern <- ifelse(status == 1, 1, ifelse(time < 0.5, 2, 3))
  delta <- rep(1, n)

  data.frame(time, status, x = rep(1, n), pattern, delta, t.time, c.time)

}

# simu_mnar_1grp(100, 1, 3.25, 1, 1, 2, 2)
# simu_mnar_1grp(100, 1, 3.25, 1, 1, 2, 2, lambda_control = 1, beta_control = 2)

simu_mnar_1grp <- function(n, delta_mnar, trunc_time, lambda, beta, lambda_cen, beta_cen, lambda_control = NULL, beta_control = NULL, seed = NULL){

  if(! is.null(seed)){
    set.seed(seed)
  }

  x <- rnorm(n,0,1)
  x <- x - mean(x)
  ## generate event time before censoring
  v  <- runif(n=n)
  t.time  <- (- log(1-v) / (lambda * exp(x * beta)))

  ## generate censoring time
  v  <- runif(n=n)
  c.time0  <- (- log(1-v) / (lambda_cen * exp(x * beta_cen)))
  c.time <- pmin(trunc_time, c.time0)

  time <- pmin(t.time, c.time)

  status <- as.numeric(t.time < c.time)

  pattern <- ifelse(status == 1, 1,
                    ifelse(time >= trunc_time, 2, 3) )

  delta <- c(1,1, delta_mnar)[pattern]

  ## generate event time after censoring
  v  <- runif(n=n)
  if(is.null(beta_control)){
    .time <- (- log(1-v) / (delta_mnar * lambda * exp(x * beta)))
  }else{
    .time <- (- log(1-v) / (delta_mnar * lambda_control * exp(x * beta_control)))
  }
  t.time <- ifelse(pattern %in% c(1,2), t.time, c.time0 + .time)

  data.frame(time, status, x, pattern, delta, t.time, c.time)

}


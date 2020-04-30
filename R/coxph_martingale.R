#' Martingale Details of a Cox Model Fit
#'
#' @inheritParams coxph_mi
coxph_martingale <- function(fit, .id = NULL, .surv = NULL, .x = NULL, .time_grid = NULL){

  if(is.null(.id)){
    .id <- fit$id
  }

  if(is.null(.surv)){
    .surv <- fit$y
  }else{
    stopifnot(class(.surv) == "Surv")
  }

  if(is.null(.x)){
    .x <- fit$x
  }else{
    .x <- as.matrix(.x)
  }

  if(is.null(.time_grid)){
    .time_grid <- fit$y[,1]
  }else{
    .time_grid <- as.numeric(.time_grid)
  }
  .time_grid <- sort(unique( c(.time_grid, .surv[,1]) ) )

  if(ncol(.x) == 1 & length(unique(.x[,1])) == 1){
    fit$coefficients <- 0
  }


  stopifnot( nrow(.x) == length(.surv))
  stopifnot( all(.surv[,1] > 0) )
  stopifnot( all(.time_grid > 0) )

  fit_db   <- cbind(fit$y[,1], fit$y[,2], fit$x)
  pred_db  <- cbind(.surv[,1], .surv[,2], .x)

  if(all(fit$id %in% .id)){
    sub <- .id %in% fit$id
  }else{
    stop("All data in fit object should be included in .surv and .x")
  }

  # Dimension
  n_fit <- nrow(fit$x)
  n   <- nrow(.x)
  n_p <- ncol(.x)
  n_t <- length(.time_grid)

  # Coxph Details (Fisher Information Matrix)
  fit_detail <- coxph.detail(fit)

  # Objects in event time scale
  fit_t_0  <- data.frame(hazard = 0, time = 0)
  fit_t    <- basehaz(fit, centered = FALSE)
  fit_t    <- rbind(fit_t_0, fit_t)

  index_to_grid <- unlist(lapply(.time_grid, function(x) sum(x >= fit_t$time)))
  fit_t_grid <- fit_t[index_to_grid, ]
  fit_t_grid$time <- .time_grid
  fit_t_grid$h0 <-  with(fit_t_grid, c(hazard[1], diff(hazard)) ) # baseline hazard

  # Objectis in subject scale
  fit_s         <- as.data.frame(as.matrix(.surv))
  fit_s$risk   <- as.numeric(exp( .x %*% fit$coefficients))

  # hazard matrix (subject by event time)
  st_hazard <- outer(fit_s$risk,fit_t_grid$h0, "*")

  # at risk process (subject by event time)
  st_y      <- outer(fit_s$time, fit_t_grid$time, ">=")
  st_dy     <- outer(fit_s$time, fit_t_grid$time, "==")

  # Conting Process matrix (subject by event time)
  st_dn  <- st_dy * fit_s$status

  # Martinal matrix (subject by event time)
  st_dm  <- st_dn * st_y - st_hazard * st_y

  # Score function matrix (subject by number of covariates)

  t_s0 <- colMeans( (fit_s$risk * st_y)[sub, ] )               # S0

  sp_score <- apply(.x, 2, function(x){
    t_s1  <- colMeans( (fit_s$risk * st_y * x)[sub, ])         # S1
    t_e   <- ifelse(t_s0 == 0, 0, t_s1 / t_s0)                 # E = S1 / S0
    score <- rowSums(outer(x, t_e, "-") * st_dm, na.rm = TRUE) # Score function
    score
  })

  # Inverse of Information Matrix (covariance p by p)
  if(n_p == 1){ imat =   sum(fit_detail$imat) / n_fit}
  if(n_p > 1){  imat = apply(fit_detail$imat, c(1,2), sum) / n_fit}

  if(ncol(.x) == 1 & length(unique(.x[,1])) == 1){
    inv_imat <- 0
  }else{
    inv_imat <- solve(imat)   # A^-1
  }


  list(sub   = sub,
       fit_s = fit_s,
       fit_t = fit_t_grid,
       st_hazard = st_hazard,
       st_dm = st_dm,
       t_s0 = t_s0,
       sp_score = sp_score,
       inv_imat = inv_imat)
}

#' Martingale Details of a Cox Model Fit
#'
#' @inheritParams coxph_mi
coxph_martingale_validate <- function(fit){


  n   <- nrow(fit$x)
  n_p <- ncol(fit$x)
  n_t <- length(unique(fit$y[,1]))

  fit_detail <- coxph.detail(fit)

  # Objects in event time scale
  fit_t    <- basehaz(fit, centered = FALSE)
  fit_t$h0 <-  with(fit_t, c(hazard[1], diff(hazard)) ) # baseline hazard

  # Objectis in subject scale
  fit_s         <- as.data.frame(as.matrix(fit$y))
  fit_s$risk   <- predict(fit, type = "risk")

  # hazard matrix (subject by event time)
  st_hazard <- outer(fit_s$risk,fit_t$h0, "*")

  # at risk process (subject by event time)
  st_y      <- outer(fit_s$time, fit_t$time, ">=")
  st_dy     <- outer(fit_s$time, fit_t$time, "==")

  # Conting Process matrix (subject by event time)
  st_dn  <- st_dy * fit_s$status

  # Martinal matrix (subject by event time)
  st_dm  <- st_dn * st_y - st_hazard * st_y

  # Score function matrix (subject by number of covariates)

  t_s0 <- colMeans(fit_s$risk * st_y)                            # S0

  sp_score <- apply(fit$x, 2, function(x){
    t_s1  <- colMeans(fit_s$risk * st_y * x)            # S1
    t_e   <- t_s1 / t_s0                                # E = S1 / S0
    score <- rowSums(outer(x, t_e, "-") * st_dm)        # Score function
    score
  })

  # Inverse of Information Matrix (covariance p by p)
  if(n_p == 1){ imat =   sum(fit_detail$imat) / n}
  if(n_p > 1){  imat = apply(fit_detail$imat, c(1,2), sum) / n}

  inv_imat <- solve(imat)   # A^-1

  list(fit_s = fit_s,
       fit_t = fit_t,
       st_hazard = st_hazard,
       st_y = st_y,
       st_dy = st_dy,
       st_dm = st_dm,
       t_s0 = t_s0,
       sp_score = sp_score,
       inv_imat = inv_imat)
}

# Validation Mode
seed <- 123
validate <- TRUE
#######################

par0.delta<-par.delta<- 2
lambda    <- 0.35
beta      <- c(0.75)
lambdaC   <- 0.15
betaC     <- c(0.75)
LL        <- 3
n.eval    <- 1000
n.mi      <- 3
n.wb      <- 100
n=500

set.seed(seed)

tau<-3.25#runif(n,1.5,3)
eval.times<-seq(0,LL,length.out = n.eval)
# true<-fctn_true(n=10000,lambda,beta,lambdaC,betaC,par0.delta,eval.times,LL,n.eval) ## true survival curve
# tauLL<-sum(true*(LL/n.eval)) ## true restricted mean survival time
#
## generate time-indept covariate
dimX<-1
X<-rnorm(n,0,1)
X<-matrix(X,n,dimX)
xx1<-matrix(X)
xx1mean<-apply(xx1,2,mean)
xx1[,1]<-xx1[,1]-xx1mean[1]
X<-xx1

## generate T
v  <- runif(n=n)
T  <- (- log(1-v) / (lambda * exp(X %*% beta)))
## generate C
v  <- runif(n=n)
C  <- (- log(1-v) / (lambdaC * exp(X %*% betaC)))
DeltaT <- (T<tau ) & (T<C)
#summary(T)
#summary(C)
#mean(DeltaT)
R <- (T>=tau)*(C>=tau) +((T>=C)*(C<tau))*2
U<-DeltaT*T+(R==1)*tau+(R==2)*C


# fctn_est<-function(U,R,n,xx1,eval.times,par.delta){
#summary(U)
deltaT <-1*(R==0)
deltaC <-1*(R>0)
# table(R)
# table(R)/n
## observed treatment
data1      <- list(time=U,status=deltaT,xx1=xx1)
fit <- coxph(Surv(time, status) ~ xx1, data1, x = TRUE, y = TRUE)


f1 <- coxph_martingale(fit)
f2 <- coxph_martingale_validate(fit)

test_that("Validate coxph_martingale without new data",{
  for(i in 1:length(f1)){
    expect_equivalent(f1[[i]], f2[[i]])
  }
})

# library(survival)
# library(testthat)

validate <- TRUE
seed <- 1234

par.delta<- 1
tau <- 3.25
LL <- 3
n.eval <- 1000

# Data Simulation
n <- 500
t.time <- rexp(n)   # event time
c.time <- runif(n, min = 0.1, max = 5)  # censoring time
c.time <- pmin(c.time, tau)

time <- pmin(t.time, c.time)
status <- as.numeric(t.time < c.time)
pattern <- ifelse(status == 1, 1, ifelse(time < 0.5, 2, 3))


# KM estimator

fit_km <- survfit(Surv(time, status) ~ 1)

fit <- coxph(Surv(time, status) ~ rep(1, n))

expect_equivalent(fit_km$surv, survfit(fit)$surv, tolerance = 3)

plot(survfit(fit)$time, survfit(fit)$surv, type = "l")
lines(fit_km$time, fit_km$surv, col = 2, type = "l")

# Match Parameters
n.mi <- 3
n.wb <- 5

## generate T
T  <- t.time


eval.times<-seq(0,LL,length.out = n.eval)
## generate C
C  <- c.time
DeltaT <- (T<tau ) & (T<C)
#summary(T)
#summary(C)
#mean(DeltaT)
R <- (T>=tau)*(C>=tau) +((T>=C)*(C<tau))*2
U<-DeltaT*T+(R==1)*tau+(R==2)*C

deltaT <-1*(R==0)
deltaC <-1*(R>0)

betahat1   <- 0
xx1 <- 1
#v.betahat1 <- vcov(fit)
#betahat1

# Main code

#fit$linear.predictors
ss<-survfit(fit)#Compute the survival function
cumu1.hazard<- -log(ss$surv)#base cummulative hazard
obs1.times  <- ss$time
base.hazard0<- cumu1.hazard-c(0,cumu1.hazard[1:(length(cumu1.hazard)-1)])#baseline hazard??
hazard0     <- base.hazard0
ltime       <- length(obs1.times)

#length(hazard0)
matobs.times<- matrix(obs1.times,nrow=n,ncol=ltime,byrow=TRUE)
matxx1      <- matrix(xx1,nrow=n,ncol=ltime,byrow=FALSE)
matdeltaT   <- matrix(deltaT,nrow=n,ncol=ltime,byrow=FALSE)
matdeltaC   <- matrix(deltaC,nrow=n,ncol=ltime,byrow=FALSE)
matU        <- matrix(U,nrow=n,ncol=ltime,byrow=FALSE)
mathazard0  <- matrix(hazard0,nrow=n,ncol=length(hazard0),byrow=TRUE)
matR        <- matrix(R,nrow=n,ncol=ltime,byrow=FALSE)
matYu       <- matobs.times <= matU
dNu         <- (matU==matobs.times)*matdeltaT

## exp(-integrated hazard) is approx by prod of (1-hazard);

hazard1 <- piTtilde <- mathazard0*exp(matxx1*betahat1[1])
hazard1[hazard1>1] <- 1
temp1   <- 1-piTtilde
temp1[temp1<0] <- 1
S1u     <- t(apply(temp1,1,cumprod))
dA1u <- piTtilde*matYu
dMu  <- (dNu-dA1u)*matYu

# dim(dMu)
# dim(matxx1)
# sum(dMu*matYu)
# sum(matxx1*dMu*matYu)

matSx0 <-          exp(matxx1*betahat1[1])*matYu
matSx1 <- matxx1  *exp(matxx1*betahat1[1])*matYu
matSx2 <- matxx1^2*exp(matxx1*betahat1[1])*matYu
sx0    <- apply(matSx0,2,mean)
sx1    <- apply(matSx1,2,mean)
sx2    <- apply(matSx2,2,mean)
matsx0 <- matrix(sx0,nrow=n,ncol=ltime,byrow=TRUE)
matsx1 <- matrix(sx1,nrow=n,ncol=ltime,byrow=TRUE)
matsx2 <- matrix(sx2,nrow=n,ncol=ltime,byrow=TRUE)

## max partial likelihood estimator of beta
matE   <- ifelse( matsx0 != 0, matsx1/matsx0, 0) # modified to use 0 for 0/0 at tails

H1i    <- apply((matxx1-matE)*dMu,1,sum)
tempu  <- ifelse( sx0 != 0, (sx2/(sx0)-(sx1/sx0)^2), 0) # modified to use 0 for 0/0 at tails

Abeta  <-mean(apply( matrix(tempu,nrow=n,ncol=ltime,byrow=TRUE)*dNu,1,sum ))
VH <- mean(apply((matxx1-matE)^2*mathazard0*exp(matxx1*betahat1[1])*matYu,1,sum))
ve.beta1<-  Abeta^(-2)*var(H1i)/n # check for linearization of betahat
# ve.beta1

matH1i <- matrix(H1i,nrow=n,ncol=ltime,byrow=FALSE)
hazard2 <- piTtilde2 <- (matR<2)*mathazard0*exp(matxx1*betahat1[1])+
  (matR==2)* mathazard0*exp(matxx1*betahat1[1])*matYu +
  (matR==2)* mathazard0*par.delta*exp(matxx1*betahat1[1])*(1-matYu)
hazard2[hazard2>1] <- 1
temp2 <- 1-piTtilde2
temp2[temp2<0] <- 1
S2u <- t(apply(temp2,1,cumprod))

## conditional survival function
temp<-apply( S1u*(matobs.times==matU), 1,sum)
temp[which(temp==0)]<-1
S1C<-matrix( temp ,n,ltime,byrow=FALSE )
S2Cu<-((S1u/S1C)*(matR==1)+(S1u/S1C)^(par.delta)*(matR==2))*(1-matYu)

## components for the martingale series
piTtilde3 <- (matR==1)* mathazard0*exp(matxx1*betahat1[1])*matxx1*(1-matYu) +
  (matR==2)* mathazard0*par.delta*exp(matxx1*betahat1[1])*matxx1*(1-matYu)
temp3 <- 1-piTtilde3
temp3[temp3<0] <- 1
S3u <- t(apply(temp3,1,cumprod))

piTtilde4 <-(matR==1)* mathazard0*exp(matxx1*betahat1[1])*matE*(1-matYu) +
  (matR==2)* mathazard0*par.delta*exp(matxx1*betahat1[1])*matE*(1-matYu)
temp4 <- 1-piTtilde4
temp4[temp4<0] <- 1
S4u <- t(apply(temp4,1,cumprod))

forg0u      <- (1-matYu)*S2Cu*exp(matxx1*betahat1[1])
g0u         <- apply(forg0u,2,mean)
forg1u      <- -log(S3u)
g1u         <- apply((1-matYu)*S2Cu*forg1u,2,mean)
forg2u      <- -log(S4u)
g2u         <- apply((1-matYu)*S2Cu*forg2u,2,mean)

matg0u      <- matrix(g0u,nrow=n,ncol=ltime,byrow=TRUE)
matg1u      <- matrix(g1u,nrow=n,ncol=ltime,byrow=TRUE)
matg2u      <- matrix(g2u,nrow=n,ncol=ltime,byrow=TRUE)
# phi1u_left  <- (matg2u-matg1u)*(Abeta)^(-1)*matH1i
phi1u_left <- 0  # Should be 0 for MAR

matg0_sx0 <- ifelse(matsx0 != 0, matg0u/matsx0, 0)  # modified to use 0 for 0/0 at tails
phi1u_righ  <- t(apply(matg0_sx0*dMu*(1-matYu) ,1,cumsum))
phi1u       <- phi1u_left-phi1u_righ

# phi1u is 0 ???
sum(abs(phi1u))

## multiple imputation


S1u.imp <- matrix(NA,nrow=n,ncol=ltime)
matR2   <- (matR==2)
S1u.imp[(matobs.times>=matU)&(matR>0)] <- S2u[(matobs.times>=matU)&(matR>0)]
options(warn=-1)
forimpu <- apply(S1u.imp,1,max,na.rm=TRUE) # the upper bound for u's
options(warn=0)

Imp.mi  <- Imp.surv <- matrix(NA,nrow=n.mi,ncol= sum(forimpu>0))

for(mm in 1:n.mi){
  jj<-0
  for(kk in 1:n){
    if(forimpu[kk]>0){
      jj<-jj+1

      # Ensure random number is the same in validaton mode
      if(! is.null(seed)){ set.seed(seed + mm * 10000 + kk) }

      if(validate){
        ui <- forimpu[kk] / mm
      }else{
        ui<-runif(1,0,forimpu[kk])
      }

      temp<-S2u[kk,]
      temmp<-which(temp>=ui)
      if(length(temmp)>0 ){ Imp.mi[mm,jj] <-obs1.times[max(temmp)]}
      if(length(temmp)==0){ Imp.mi[mm,jj] <-obs1.times[ltime]     }
    }
  }
}

forV1.imp  <- array(NA, c(n.mi,n,ltime))
for(kk in 1:n.mi){
  forV1.imp[kk,,] <- (matrix(c(Imp.mi[kk,], U[R==0]-U[R==0]),nrow=n,ncol=ltime,byrow=FALSE)>=matobs.times)*(1-matYu)
}

## est for each imputed dataset
S1u.adj <- mi.adj <- apply(matrix(c(Imp.mi[1,], U[R==0]),nrow=n,ncol=ltime,byrow=FALSE)>=matobs.times,2,mean )
## ve for each imputed dataset
vemi.adj<- apply(matrix(c(Imp.mi[1,], U[R==0]),nrow=n,ncol=ltime,byrow=FALSE)>=matobs.times,2,var )/n #mi.adj *(1-mi.adj)/(n-1)
thisimp <- c(Imp.mi[1,], U[R==0])
thisimp[which(thisimp>LL)]<-LL
mi.tau     <- mean(thisimp)
vemi.tau   <- var(thisimp)/n

for(mm in 2:n.mi){
  S1u.adj<-S1u.adj+
    apply(matrix(c(Imp.mi[mm,], U[R==0]),nrow=n,ncol=ltime,byrow=FALSE)>=matobs.times,2,mean )
  thismi.adj<-apply(matrix(c(Imp.mi[mm,], U[R==0]),nrow=n,ncol=ltime,byrow=FALSE)>=matobs.times,2,mean )
  mi.adj  <-rbind(mi.adj, thismi.adj )
  vemi.adj<-rbind(vemi.adj, apply(matrix(c(Imp.mi[mm,], U[R==0]),nrow=n,ncol=ltime,byrow=FALSE)>=matobs.times,2,var )/n)
  thisimp<-c(Imp.mi[mm,], U[R==0])
  thisimp[which(thisimp>LL)]<-LL
  mi.tau  <- c(mi.tau, mean(thisimp) )
  vemi.tau<- c(vemi.tau, var(thisimp)/n  )

}
S1u.adj <- S1u.adj/n.mi
loc.rmst<-which( obs1.times<LL )
dur.rmst<-c(obs1.times[loc.rmst],LL)-c(0,obs1.times[loc.rmst])
tau.adj <-sum(c(1,S1u.adj[loc.rmst]) * dur.rmst)
## rubin's rule for combining
ve.adjmi<- (n.mi+1)/(n.mi)*apply(mi.adj,2,var) + apply(vemi.adj,2,mean)
ve.taumi<- (n.mi+1)/(n.mi)*var(mi.tau) + mean(vemi.tau)

## wild-bootstrap
forV1.mean <- apply(forV1.imp,c(2,3),mean)
est.wb1<-matrix(0,n.wb,ltime)
tau.wb1<-rep(0,n.wb)
for( bb in 1:n.wb ){
  for( mm in 1:n.mi){

    # Ensure random number is the same in validaton mode
    if(! is.null(seed)){set.seed(seed + bb * 10000 + mm)}

    u.wb        <- matrix( rnorm(n), nrow=n,ncol=ltime,byrow=FALSE)
    thismmu     <- apply( (forV1.imp[mm,,] - forV1.mean)*(1-matYu)*u.wb,2,sum )/n/n.mi
    est.wb1[bb,]<- est.wb1[bb,]+ thismmu
    tau.wb1[bb] <- tau.wb1[bb]+sum(c(0,thismmu[loc.rmst]) * dur.rmst)
  }
}
forV1<-apply(est.wb1,2,var)
forV1tau<- var(tau.wb1)

## smooth poi estimates
est<-approx(obs1.times, S1u.adj, eval.times, method="constant" )$y
#plot(obs1.times,S1u.adj)
#lines(eval.times,est,col=2)

## smooth var estimates
matS1u.adj <- matrix(S1u.adj,nrow=n,ncol=ltime,byrow=TRUE)
forV2      <- apply( (phi1u+matYu*(matU>=matobs.times)+(1-matYu)*S2Cu - matS1u.adj)^2 , 2 , mean )/n
forV2.wb   <-        (phi1u+matYu*(matU>=matobs.times)+(1-matYu)*S2Cu - matS1u.adj)

est.wb2<-matrix(NA,n.wb,ltime)
tau.wb2<-rep(0,n.wb)
for( bb in 1:n.wb ){

  # Ensure random number is the same in validaton mode
  if(! is.null(seed)){set.seed(seed + bb)}

  u.wb        <- matrix( rnorm(n), nrow=n,ncol=ltime,byrow=FALSE)
  est.wb2[bb,]<- apply(forV2.wb*u.wb,2,mean)
  tau.wb2[bb] <- sum(c(0,est.wb2[bb,loc.rmst]) * dur.rmst)
}
est.wb <- est.wb2+est.wb1
tau.wb <- tau.wb2+tau.wb1

ve.wb1 <- apply(est.wb,2,var)
ve.wb2 <- (apply(est.wb,2,quantile,0.975,type=5) - apply(est.wb,2,quantile,0.025,type=5) )^2/16
est.q1 <- S1u.adj-apply(est.wb,2,quantile,0.975,type=5)
est.q2 <- S1u.adj-apply(est.wb,2,quantile,0.025,type=5)
tau.q1 <- tau.adj-quantile(tau.wb,0.975,type=5)
tau.q2 <- tau.adj-quantile(tau.wb,0.025,type=5)
#spl.out2 <- smooth.spline(obs1.times, ve.adj, df = 10)
#ve <- predict( spl.out2, eval.times)$y

ve.adj <- forV2 +forV1 #ve.adj <- ve.wb1
ve<-approx(obs1.times, ve.adj, eval.times, method="constant" )$y
est.q1<-approx(obs1.times, est.q1, eval.times, method="constant" )$y
est.q2<-approx(obs1.times, est.q2, eval.times, method="constant" )$y
ve.adjmi<-approx(obs1.times, ve.adjmi, eval.times, method="constant" )$y

#var(tau.wb1)+var(tau.wb2)
ve.tau<-var(tau.wb)

# }

##################Split Wild Bootstrap###################################################


fit$id <- 1:length(time)
x <- matrix(1, ncol = 1, nrow = n)
wb_val <- coxph_wb_utility_simple(fit = fit,
                                  id = fit$id,
                                  time = as.numeric(time),
                                  status = as.numeric(status),
                                  x = x,
                                  pattern = pattern,
                                  delta = c(1,1,1)[pattern])

expect_equivalent(wb_val$phi, phi1u)
expect_equivalent(wb_val$st_delta_survival, S2u)
expect_equivalent(wb_val$st_delta_con_survival, S2Cu)



u_time <- sort(unique(time))
mi_t <- mi_time(time,
                status,
                u_time ,
                wb_val$st_delta_survival,
                n_mi = n.mi,
                seed = seed,
                validate = validate)

mi_surv <- mi_survival(time, u_time, mi_t)
mi_est_rmst <- mi_rmst(mi_t, tau = LL)
wb_var <- wild_variance(time = time,
                        status = status,
                        u_time = u_time,
                        imp_time = mi_t,
                        s_mi = mi_surv[,1],
                        phi = wb_val$phi,
                        phi_id = 1:n,
                        st_con_survival = wb_val$st_delta_con_survival,
                        n_b = n.wb,
                        tau = LL,
                        seed = seed,
                        validate = TRUE)

test_that("Validate Multiple Imputation", {
  expect_equivalent(mi_surv[,1], S1u.adj)
  expect_equivalent(mi_surv[,2], sqrt((n.mi+1)/(n.mi)*apply(mi.adj,2,var) + apply(vemi.adj,2,mean)))
  expect_equivalent(mi_est_rmst[2], sqrt(ve.taumi))
})

test_that("Validate Wild Bootstrap", {
  expect_equivalent(wb_var$surv_wb_sd, sqrt(ve.adj))
  expect_equivalent(wb_var$rmst_wb_sd, sqrt(ve.tau))
})


tmp <- rmst_delta(time = time,
                  status = status,
                  x = x, group = rep(1, n), pattern, c(1,1,1)[pattern],
                  tau = LL, n_mi = n.mi, n_b = n.wb, seed = seed, wild_boot = TRUE, validate = TRUE)


mi_rmst_val <- c(tau.adj, sqrt(ve.taumi), sqrt(ve.tau))
mi_s_val <- cbind(s = S1u.adj,
                  mi_sd = sqrt((n.mi+1)/(n.mi)*apply(mi.adj,2,var) + apply(vemi.adj,2,mean)),
                  wb_sd = sqrt(ve.adj)
)

# expect_equivalent(tmp$fit_wb$st_delta_survival, wb_val$st_delta_survival)
# expect_equivalent(tmp$fit_wb$st_delta_con_survival, wb_val$st_delta_con_survival)
# expect_equivalent(tmp$fit_wb$phi, wb_val$phi)
# expect_equivalent(tmp$wb_var$surv_wb_sd, wb_var$surv_wb_sd)

expect_equivalent(tmp$rmst[, c("sd", "wb_sd")], mi_rmst_val[-1])
expect_equivalent(tmp$surv[[1]][,c(-1, -2)], mi_s_val)


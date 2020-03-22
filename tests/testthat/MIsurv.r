rm(list=ls())
library(zoo)
library(MASS)
library(survival)
library(parallel)
library(rootSolve)
library(plotrix)
source("fctn_MIsurv.r")
set.seed(1234)
main
numWorkers <- 16
ranseed <- 1:2000
sim<-length(ranseed)
mc=1

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

pool.delta<-c(1,1.5,2,2.5,3)
n.col<-length(pool.delta)
cname<-c("true","est","sd","se.mi","reb.mi","cvg.mi","se.wild","reb.wild","cvg.wild")
n.row<-length(cname)
Res.table <-matrix(NA,nrow=n.row,ncol=n.col)
colnames(Res.table)<-pool.delta
row.names(Res.table)<-cname

char1<-paste(paste(paste(paste(paste("MIsurv.","n",sep=""),n,sep=""),".mi",sep=""),n.mi,sep=""),".pdf",sep="")

pdf(char1)
jj<-main(10001)
plot(jj$eval.times,jj$true,ylim=c(0.2,1),type="n",ylab="S(t)",xlab="t")


#############################################################################################################
for(kk in 1:n.col){
  
  par.delta <- pool.delta[kk]
  
  if(mc==1){
    res <- mclapply(ranseed,main,mc.cores=numWorkers)
  }
  if(mc==0){
    res <- lapply(ranseed,main)
  } 
  
  char<-paste(paste(paste(paste(paste(paste(paste("MIsurv.","n",sep=""),n,sep=""),".mi",sep=""),
                                n.mi,sep=""),".delta",sep=""),par.delta,sep=""),".Rdata",sep="")
  char
  #load(char)
  save.image(char)
  
  eval.times<- res[[1]]$eval.times
  betahat1  <- res[[1]]$betahat1
  ve.beta1  <- res[[1]]$ve.beta1
  
  est1<- res[[1]]$est1
  se1.mi<- res[[1]]$se1.mi
  se1 <- res[[1]]$se1
  true      <- res[[1]]$true
  est.q1    <- res[[1]]$est.q1
  est.q2    <- res[[1]]$est.q2
  
  est2<- res[[1]]$est2
  se2.mi <- res[[1]]$se2.mi
  se2 <- res[[1]]$se2
  tauLL     <- res[[1]]$tauLL 
  tau.q1    <- res[[1]]$tau.q1
  tau.q2    <- res[[1]]$tau.q2
  
  cnt<-0
  for(i in 2:sim){
    if(class(res[[i]]) == "try-error"){cnt<-cnt+1;next}
    betahat1  <- rbind(betahat1,res[[i]]$betahat1)
    ve.beta1  <- rbind(ve.beta1,res[[i]]$ve.beta1)
    est1<- rbind(est1,res[[i]]$est1)
    se1 <- rbind(se1,res[[i]]$se1)
    se1.mi<- rbind(se1.mi,res[[i]]$se1.mi)
    true      <- rbind(true,res[[i]]$true)
    est.q1    <- rbind(est.q1,res[[i]]$est.q1)
    est.q2    <- rbind(est.q2,res[[i]]$est.q2)
    
    est2<- rbind(est2,res[[i]]$est2)
    se2 <- rbind(se2,res[[i]]$se2)
    se2.mi <- rbind(se2.mi,res[[i]]$se2.mi)
    tauLL     <- rbind(tauLL,res[[i]]$tauLL)
    tau.q1    <- rbind(tau.q1,res[[i]]$tau.q1)
    tau.q2    <- rbind(tau.q2,res[[i]]$tau.q2)
  }
  nsim<-dim(est1)[1]
  
  dim(est1)
  apply(betahat1,2,mean)
  apply(betahat1,2,var  )
  apply(ve.beta1,2,mean)
  
  # summary 1
  estmean <- apply(est1,2,mean)
  estse   <- apply(se1,2,mean)
  trueval <- apply(true,2,mean)
  truemat <- matrix( trueval,n=nsim,ncol=length(trueval),byrow=TRUE)
  estbias <- apply(est1-truemat,2,mean )
  left    <- (est1-2*(se1))<truemat
  right   <- (est1+2*(se1))>truemat
  cvg     <- left*right
  
  left.mi <- (est1-2*(se1.mi))<truemat
  right.mi<- (est1+2*(se1.mi))>truemat
  cvg.mi  <- left.mi*right.mi
  
  left1   <- (est.q1)<truemat
  right1  <- (est.q2)>truemat
  cvg1    <- left1*right1
  
  # bias
  round( (estmean-trueval)*100,0)
  # sd
  sd1<-apply(est1,2,sd  )
  sd1hat<-apply(se1 ,2,mean)
  sd1hat.mi<-apply(se1.mi ,2,mean)
  round(apply(est1,2,sd  )*100,2)
  round(apply(se1 ,2,mean)*100,2)

  q1   <- estmean-2*sd1
  q2   <- estmean+2*sd1
  q1hat<- estmean-2*sd1hat
  q2hat<- estmean+2*sd1hat
  q1mi <- estmean-2*sd1hat.mi
  q2mi <- estmean+2*sd1hat.mi
  
  
  
  ############################
  Res.table["true"    ,kk]=mean(tauLL)
  Res.table["est"     ,kk]=mean(est2)
  Res.table["sd"      ,kk]=sd(est2)
  Res.table["se.wild" ,kk]=mean(se2)
  Res.table["se.mi"   ,kk]=mean(se2.mi)
  
  Res.table["reb.wild",kk]=round((mean(se2)/sd(est2)-1)*100,1)
  Res.table["reb.mi"  ,kk]=round((mean(se2.mi)/sd(est2)-1)*100,1)
  
  left0    <- (est2-2*(se2))<tauLL
  right0   <- (est2+2*(se2))>tauLL
  cvg0     <- left0*right0
  Res.table["cvg.wild",kk]=mean(cvg0)
  
  left0.mi    <- (est2-2*(se2.mi))<tauLL
  right0.mi   <- (est2+2*(se2.mi))>tauLL
  cvg0.mi     <- left0.mi*right0.mi
  Res.table["cvg.mi"  ,kk]=mean(cvg0.mi)
  
  ############################
  lines(eval.times,trueval)
  lines(eval.times,estmean,col=kk+1)
  
  
}

legend("top", legend=paste(expression("delta="),pool.delta,sep=""),
       col=(1:n.col)+1, lty=1)

dev.off()


Res.table


rm(list = ls())
######################### Type1 error / Power #########################
setwd('D:\\Desktop\\code')
library(ggplot2)
library(corpcor)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(MTS)
sourceCpp('Online.cpp')
sourceCpp('SN_Online.cpp')
source('gendataMulLast.R')
## boundary function 
gamma<-0
Qk<-function(s0,gamma0){
  u=(1+s0)*(s0/(1+s0))^gamma0
}


DGP <- 204 ## Data generation 


#################### para
T <- 100 ## Trian_Phase
L <- 1   ## monitoring period


# yt <- genDataMul(L,T,model = DGP)
# MTSplot(yt)
d <- 2 ## Data dimension

kernal <- "gauss" ##  gauss
delta0 <- 1

###### The weight of the statistic ###########
Tm <- T
TypeI <- 0.05
N<-(1+L)*T
k<-1:(N-T)
## boundary function
vaha<-(Tm+k)^2/(Tm*Qk(k/Tm,gamma)^2)



#######  Bootstrap ################## 

KT <- max(5, sqrt(log10(T)))  ## select m^
## flat-top lag-window
lambda <- function(t) {
  abst <- abs(t)
  (abst <= .5) + 2*(1 - abst)*((abst > .5) & (abst <= 1))
}

nsim<-1000
Tn <- matrix(0, nrow=nsim, ncol=1) ## Calculate the statistics for each generated data
Tn.star <- matrix(0, nrow=nsim, ncol=1) ## Calculate statistics for each Bootsrap data
rej0 <- rep(0,1) ## reject rate

startTime <- Sys.time()
for (s in 1:nsim) {
  
  yt <- genDataMul(L,T,model = DGP)
  Block <- numeric(2)
  for (i in 1:d) {
    ut <- yt[1:T,i]
    R <- as.vector(acf(ut, lag.max=2*KT, plot=FALSE)[[1]])
    tmp <- which(abs(R[1:KT]) < 2*sqrt(log10(T)/T)) 
    if (length(tmp) > 0) {
      M <- 2*(tmp[1] - 1)
    } else {
      M <- 2*(KT - 1)
    }
    ghat <- sum(lambda((-M:M)/M)*R[abs(-M:M)+1])
    Ghat <- sum(lambda((-M:M)/M)*R[abs(-M:M)+1]*abs(-M:M))
    D.SB <- 2*ghat^2
    blocksize <- max(1, round((2*Ghat^2/D.SB)^(1/3)*T^(1/3)))
  Block[i] <- blocksize
  }
  #Block
  blocksize <- round( mean(Block) )
  # blocksize
  
  y.star <- matrix(rep(NA,1*2),nrow=1,ncol=d)

  while (nrow(y.star) < (N+T/2)) {
    idx.start <- sample.int(T, 1, replace=TRUE) 
    idx.len <- rgeom(1, 1/blocksize) + 1       
    idx <- idx.start:(idx.start + idx.len)
    idx <- (idx - 1) %% T + 1
    idx
    y.star <- rbind(y.star, yt[(idx),])
  }
  y.star <- y.star[-1,]
  y.star <- y.star[1:N,]
  
  
  Tn[s] <- max( statMarginal(yt,T,1,vaha,kernal,delta0) )  
  Tn.star[s] <- max( statMarginal(y.star,T,1,vaha,kernal,delta0) ) 
  ## statMarginal  stands for D_m(k)
  ## SN_statMarginal  stands for D_m^{SN}(k)
  
  if (s %% 50 == 0) {
    critvals <-  quantile(Tn.star[1:s,], 1-TypeI)
    rej0<- sum(critvals < Tn[1:s])
    
    cat('\n')
    rej <- round(100*rej0/s,1)
    print(rej)
  }
  elapsed = difftime(Sys.time(),startTime,units='mins')
  remaining = round((nsim-s)*elapsed/s)
  cat('\r', floor(s/nsim*100),'% complete (',remaining,' min remaining)    ',sep='')
  flush.console()
}



critvals

rej0
rej


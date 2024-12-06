genDataMul <- function(L,T,model) {
  library(mvtnorm)
  library(MTS)
  library(corpcor)
  ### 
  burnin <- ceiling(T/2)
  N <- (1+L)*T
  obs <- burnin + N
  #cpt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
  ## parameter
  delta <- 0
  kor0 <- 0
  kor1<-0.2
  kor2<-0.6
  kor3<-0.9
  mu0=rep(0,2)
  mu1=rep(1,2)
  ##  covariance matrix
  sigma0 <- matrix(c(1,kor0,kor0,1),nrow<-2)
  sigma1 <- matrix(c(1,kor1,kor1,1),nrow<-2)
  sigma2 <- matrix(c(1,kor2,kor2,1),nrow<-2)
  sigma3 <- matrix(c(1,kor3,kor3,1),nrow<-2)
  ## VAR(1) coefficients
  B <- matrix(c(0.2,0.1,0.1,0.2),2,2)/2
  #eigen(B)$values
  ## VAR(1) coefficients 
  A1 <- matrix(c(0.5,0.2,0.2,0.1),2,2)
  #eigen(A1)$values
  ## BEKK-GARCH(1,1)
  C <- t(matrix(c(0.004,0.005,0,0.003),nrow<-2))
  Amat <- matrix(c(0.254,0.004,-0.004,0.332),nrow<-2)
  Bmat <- matrix(c(0.941,-0.019,0.023,0.864),nrow<-2)
  intercept <- t(C)%*%(C)
  
  ## If the absolute values of the matrix eigenvalues are all less than 1, then the VAR (1) model is stationary.
  reps0 <-  function(n,kr=kor) {rmvnorm(n, mean = rep(0,2), sigma = matrix(c(1,kr,kr,1),nrow=2))}
  
  reps1 <-  function(n,kr=kor) {rmvt(n, sigma = matrix(c(1,kr,kr,1),nrow = 2), df = 5)} ## for model 304
  
  reps2 <-  function(n,kr=kor) {rmvt(n, sigma = (3/5)*matrix(c(1,kr,kr,1),nrow = 2), df = 5)} #t5(0,sigma)
  
  x <- matrix(0,obs,2)
  epsN <- reps0(obs,kor0)
  epsT <- reps1(obs,kor2) ## for model 304
  
  if (model == 1) {           ##### DGP N1 #######
    x <- epsN
  } 
  else if (model == 2) {      ##### DGP N2 #######  VAR(1)
    x<-ts(VARMAsim(obs,arlags=c(1),phi=0.3*diag(2),sigma=sigma0)$series)
  }
  else if (model == 3) {      ##### DGP N3 #######  VAR(1)
    x<-ts(VARMAsim(obs,arlags=c(1),phi=0.5*diag(2),sigma=sigma0)$series)
  }
  else if (model == 4) {      ##### DGP N4 #######  VAR(1)
    B <- matrix(c(0.2,0.1,0.1,0.2),2,2)/2
    #eigen(B)$values
    x<-ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=diag(2) )$series)
  }
  else if (model == 5) {      ##### DGP N5 #######  VAR(1)
    x<-ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=sigma1)$series)
  }
  else if (model == 6) {      ##### DGP N6 #######  BEKK-GARCH(1,1)
    eps1 <- reps0(obs,kor0)
    x <- matrix(0,obs,2)
    H <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    H[[1]] <- 0.01*diag(2)      # initial 
    x[1,] <- t(chol(H[[1]]))%*%eps1[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( x[(i-1),])
      Inn <- Ri%*%t(Ri)  
      H[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%H[[(i-1)]]%*%Bmat
       h12 <-  t(chol(H[[i]]))    
      x[i,] <- h12%*%eps1[i,]
    }
  }
  
  ######################## mean change ################################
  else if (model == 101) {    #### DGP M1 δ=1 ####### 
    delta <- 1
    x <- epsN
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
    
  } 
  else if (model == 102) {     ##### DGP M2 δ=1 ####### 
    delta <- 1
    x<-ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=diag(2) )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
  } 
  else if (model == 103) {     ##### DGP M3 δ=1 ####### 
    delta <- 1
    x <- ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=diag(2))$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <- x[change.pt:obs,] + delta*mu1
  } 
  
  
  ######################## second-order moment change ################################ 
  
  ########################## Change in  Variance
  else if (model == 201) {   ##### DGP V1 δ=2    #######   
    delta <- 2
    x <- rmvnorm(obs, mean = rep(0,2), sigma = diag(2))
    y <-rmvnorm(obs, mean = rep(0,2), sigma = (1+delta)*diag(2))
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  else if (model == 202) {   ##### DGP V2 δ=2     #######   Variance
    delta <- 2
    x <- ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=diag(2) )$series)
    y <-  ts(VARMAsim(obs,arlags=c(1),phi=B,sigma=(1+delta)*diag(2) )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  else if (model == 203) {   ##### DGP V3 δ=2     #######   Variance
    delta <- 2
    x <- ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=diag(2) )$series)
    y <-  ts(VARMAsim(obs,arlags=c(1),phi=A1,sigma=(1+delta)*diag(2) )$series)
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  else if (model == 204) {   ##### DGP V5 δ=(1,1)     #######   Variance
    delta <- 1
    ### before change
    eps1 <- rmvnorm(obs, mean = rep(0,2), sigma =matrix(c(1,0,0,1),nrow<-2)  )
    
    x <- matrix(0,obs,2)
    H <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    H[[1]] <- 0.01*diag(2)      # initial 
    x[1,] <- t(chol(H[[1]]))%*%eps1[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( x[(i-1),])
      Inn <- Ri%*%t(Ri)  
      H[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%H[[(i-1)]]%*%Bmat
      h12 <-  t(chol(H[[i]]))    
      x[i,] <- h12%*%eps1[i,]
    }
    ### after change
    eps2 <-  rmvnorm(obs, mean = rep(0,2), sigma =matrix(c(1+delta,0,0,1+delta),nrow<-2)  )
    y <- matrix(0,obs,2)
    Hy <- vector("list", obs) # 创建一个空列表H,存放N个条件协方差矩阵
    Hy[[1]] <- 0.01*diag(2)      # initial 
    y[1,] <- t(chol(Hy[[1]]))%*%eps2[1,] 
    
    for (i in 2:obs) {
      Ri <- as.matrix( y[(i-1),])
      Inn <- Ri%*%t(Ri)  
      Hy[[i]] <- intercept+t(Amat)%*%Inn%*%(Amat)+t(Bmat)%*%Hy[[(i-1)]]%*%Bmat
      h12 <-  t(chol(Hy[[i]]))    
      y[i,] <- h12%*%eps2[i,]
    }
    change.pt <- floor(burnin + T*(1 + L*runif(1, 0, 4/5)))
    x[change.pt:obs,] <-y[change.pt:obs,]
  }
  x <- x[-(1:burnin),]
}


## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(knitr)

## -----------------------------------------------------------------------------
giccdmm <- function(y, z, lambda, constr=TRUE,maxit=c(20, 50)) {
  n<-length(y)
  p<-ncol(z)
  l<-length(lambda)
  gic<-rep(0,l)
  beta.total<-matrix(0,ncol=p,nrow=l)
  for(i in 1:l){
    
    if(constr){
      mu<-1
      beta.old<-rep(1,p)
      beta<-rep(0,p)
      beta0<-rep(1,p)
      alpha<-0
      lam<-lambda[i]
      k<-0
      while (k<maxit & sum(abs(beta-beta0))>10e-8) {
        beta0<-beta
        while (sum(abs(beta-beta.old))>10e-8) {
          beta.old<-beta
          for (j in 1:p) {
            t<-t(z[,j])%*%(y-z[,-j]%*%beta[-j])/n-mu*(sum(beta[-j])+alpha)
            beta[j]<-sign(t)*max(c(abs(t)-lam,0))/(sum(z[,j]^2)/n+mu)
          }
        }
        alpha<-alpha+sum(beta)
        k<-k+1
      }
      beta.total[i,]<-beta
      gic[i]<-log(sum((y-z%*%beta)^2))+(length(which(beta!=0))-1)*log(log(n))*log(max(p,n))/n
      
    }
    
    else{
      mu<-0
      beta.old<-rep(1,p)
      beta<-rep(0,p)
      beta0<-rep(1,p)
      alpha<-0
      lam<-lambda[i]
      k<-0
      while (k<maxit & sum(abs(beta-beta0))>10e-8) {
        beta0<-beta
        while (sum(abs(beta-beta.old))>10e-8) {
          beta.old<-beta
          for (j in 1:p) {
            t<-t(z[,j])%*%(y-z[,-j]%*%beta[-j])/n-mu*(sum(beta[-j])+alpha)
            beta[j]<-sign(t)*max(c(abs(t)-lam,0))/(sum(z[,j]^2)/n+mu)
          }
        }
        alpha<-alpha+sum(beta)
        k<-k+1
        
      }
      
      beta.total[i,]<-beta
      gic[i]<-log(sum((y-z%*%beta)^2))+(length(which(beta!=0))-1)*log(log(n))*log(max(p,n))/n
      
    }
    
  }
  s<-which.min(gic)
  
  return(list(beta=beta.total,gic=gic,lam=lambda))
}

## -----------------------------------------------------------------------------
library(mvtnorm)
sim.data <- function(n=50, p=30, rho=0.2, seed=123) {
	set.seed(seed)
	bet <- c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2, rep(0, p - 8))
	alp <- c(rep(p/2, 5), rep(1, p - 5))
	sig <- matrix(1, p, p)
	sig <- rho^abs(row(sig) - col(sig))
  x <- switch("ln", ln={x <- exp(rmvnorm(n, log(alp), sig)); x/rowSums(x)}, dir=rdirichlet(n, alp))
	z <- log(x)
	y <- rnorm(n, z %*% bet, 0.5)
	data <- list(y=y, z=z, bet=bet)
	
	return(data)
}
data<-sim.data()

## -----------------------------------------------------------------------------
lambda<-seq(0.05,0.15,0.02)
bet <- c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2, rep(0, 30 - 8))
result1<-giccdmm(data$y,data$z,lambda,constr = T,maxit = 40)
result2<-giccdmm(data$y,data$z,lambda,constr = F,maxit = 40)
beta1<-result1$beta[which.min(result1$gic),]
beta2<-result2$beta[which.min(result2$gic),]
PE1<-sum((data$y-data$z%*%beta1)^2)/50
PE2<-sum((data$y-data$z%*%beta2)^2)/50
loss1<-sum((bet-beta1)^2)
loss2<-sum((bet-beta2)^2)
knitr::kable(rbind(c("Type","PE","LOSS"),c("lasso with zero-sum constraint",round(c(PE1,loss1),4)),c("lasso without zero-sum constraint",round(c(PE2,loss2),4))))

## -----------------------------------------------------------------------------
subgroup_l0<-function(y,X,lam){
  n<-length(y)
  p<-ncol(X)
  B <- matrix(0, nrow = n*(n-1)/2, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      B[(2*n-i)*(i-1)/2+j-i,i]<-1
      B[(2*n-i)*(i-1)/2+j-i,j]<--1
    }
  }
  x<-solve(t(X)%*%X)%*%t(X)
  delta<-1e-5
  gamma<-2
  mu.total<-vector()
  BIC<-vector()
  beta.total<-vector()
  for (i in 1:(length(lam))) {
    lambda<-lam[i]
    beta0<-matrix(lm(y~X)$coefficient[2:(p+1)], nrow = p, ncol = 1)
    mu0<-matrix(lm(y~X)$coefficient[1], nrow = n, ncol = 1)
    w<-rep(1,n*(n-1)/2)
    f<-1
    while(f==1){
      mu<-solve(lambda*t(B)%*%diag(w)%*%B+diag(n))%*%(y-X%*%beta0)
      if(sum((mu-mu0)^2)<10-5) f<-0
      mu0<-mu
      Bmu<-B%*%mu
      for(i in 1:(n*(n-1)/2)){
        if(Bmu[i,1]<=delta) w[i]=delta^(-2)*exp(-2*log(1+abs(Bmu[i,1]/delta)^gamma)/gamma)
        else w[i]=abs(Bmu[i,1])^(-2)*exp(-2*log(1+abs(delta/Bmu[i,1])^gamma)/gamma)
      }
    }
    
    t<-table(round(mu))
    c<-which(as.numeric(t)>=0.05*n)
    mu<-round(mu)
    name<-as.numeric(names(t))
    beta.group<-matrix(0,nrow=length(which(as.numeric(t)>=0.05*n)),ncol=p)
    for(i in c){
      beta.group[which(c==i),]<-lm(y[which(mu==name[i]),]~X[which(mu==name[i]),])$coefficient[-1]
    }
    beta.group[is.na(beta.group)]<-0
    beta<-colSums( beta.group*as.numeric(t)[which(as.numeric(t)>=0.05*n)])/sum(as.numeric(t)[which(as.numeric(t)>=0.05*n)])
    
    if(sum((beta-beta0)^2)<10e-5) F<-0
    beta0<-beta
    
    F<-1
    
    while(F==1){
      
      f<-1
      k<-0
      while(f==1){
        k<-k+1
        mu<-solve(lambda*t(B)%*%diag(w)%*%B+diag(n))%*%(y-X%*%beta0)
        if(sum((mu-mu0)^2)<10-5) f<-0
        mu0<-mu
        Bmu<-B%*%mu
        for(i in 1:(n*(n-1)/2)){
          if(Bmu[i,1]<=delta) w[i]=delta^(-2)*exp(-2*log(1+abs(Bmu[i,1]/delta)^gamma)/gamma)
          else w[i]=abs(Bmu[i,1])^(-2)*exp(-2*log(1+abs(delta/Bmu[i,1])^gamma)/gamma)
        }
      }
      
      
      
      t<-table(round(mu))
      c<-which(as.numeric(t)>=0.05*n)
      mu0<-round(mu)
      name<-as.numeric(names(t))
      beta.group<-matrix(0,nrow=length(which(as.numeric(t)>=0.05*n)),ncol=p)
      for(i in c){
        beta.group[which(c==i),]<-lm(y[which(mu0==name[i]),]~X[which(mu0==name[i]),])$coefficient[-1]
      }
      beta.group[is.na(beta.group)]<-0
      beta<-colSums( beta.group*as.numeric(t)[which(as.numeric(t)>=0.05*n)])/sum(as.numeric(t)[which(as.numeric(t)>=0.05*n)])
      
      if(sum((beta-beta0)^2)<10e-5) F<-0
      beta0<-beta
    }
    
    beta.total<-rbind(beta.total,beta)
    mu.total<-cbind(mu.total,round(mu,1))
    
    BIC<-c(BIC,log(sum((y-X%*%beta-round(mu))^2)/n) +  log(log(n+2))*log(n)*sum(as.numeric(table(round(mu)))[which( as.numeric(table(round(mu)))<0.1*n)])/n)
    
    
  }
  return(list(beta=beta.total,mu=mu.total,BIC=BIC,lambda=lambda))
}
  

## -----------------------------------------------------------------------------
set.seed(12)
n<-80
library(MASS)
y<-matrix(0,ncol=1,nrow=n)

X<-mvrnorm(n,rep(0,2),diag(2))
beta<-matrix(c(3,2), nrow = 2, ncol = 1)
y<-X%*%beta+rnorm(n,0,1)

y[1:40,]<-y[1:40,]+rep(-3,40)
y[41:80]<-y[41:80]+rep(3,40)

## -----------------------------------------------------------------------------
lam<-c(0.015,0.03,0.07,0.1,0.3)
result<-subgroup_l0(y,X,lam)
result$BIC

## -----------------------------------------------------------------------------
plot(result$mu[1,],type = "l",ylim = c(-4,4),ylab = expression(mu),xaxt="n",xlab=expression(lambda))
axis(1,1:5,lam)
for (i in 2:80) {
  lines(result$mu[i,])
}


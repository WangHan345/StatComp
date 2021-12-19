## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(knitr)

## -----------------------------------------------------------------------------
data("mtcars")
set.seed(1234)
knitr::kable(mtcars[sort(sample(1:nrow(mtcars),12)),])

## -----------------------------------------------------------------------------
plot(mtcars$wt,mtcars$mpg,pch=16,xlab = "wt",ylab = "mpg")
abline(lm(mtcars$mpg~mtcars$wt),col="red")
plot(mtcars$wt,mtcars$disp,pch=16,xlab = "wt",ylab = "disp")
abline(lm(mtcars$disp~mtcars$wt),col="red")

## -----------------------------------------------------------------------------
rayleigh<-function(sigma){#inverse transform algorithm 
  n<-10000
  u<-runif(n)#Uniform(0,1)
  y<--2*log(1-u)#inverse transform of exponential distribution
  x<-sigma*sqrt(y)#generate random samples from the Rayleigh distribution
  hist(x, prob=T,main=bquote('sigma='*.(sigma)),nclass=20)
  y<-seq(0,ceiling(max(x)),0.01)
  lines(y,y*exp(-y*y/2/sigma^2)/sigma^2)
}
sigma<-c(1,3,5,10,20,30)
for (i in 1:length(sigma)) {
  rayleigh(sigma[i])
}

## -----------------------------------------------------------------------------
loc.mix<-function(p1){#generate 1000 samples from the mormal location mixture
  n1<-rbinom(1,1000,p1)
  n2<-1000-n1
  x1<-rnorm(n1,0,1)
  x2<-rnorm(n2,3,1)
  x<-c(x1,x2)
  hist(x, prob=T,main=bquote('p1='*.(p1)),nclass=20)
}
cat('p1=0.75')
loc.mix(0.75)

for (p in seq(from = 0.1, to = 0.9, by = 0.1)) {
  loc.mix(p)
}

## -----------------------------------------------------------------------------
poissongamma<-function(lambda,alpha,beta){#lambda is the parameter of Poisson process, alpha and beta are the parameters of Gamma distribution
  n<-10000#number of samples
  x<-rep(0,n)
  t<-10
  for(i in 1:n){
    y<-rgamma(1000,alpha,beta)
    T<-rexp(1000,lambda)#interval times
    S<-cumsum(T)#arrival times
    N<-min(which(S>t))-1
    x[i]<-sum(y[1:N])
    
  }
  return(c(mean(x),var(x),lambda*10*alpha/beta,lambda*10*(alpha+alpha^2)/beta^2))#output mean and variance of X(t) and lambda*t*E[Y] and lambda*t*E[Y^2]
}
p<-poissongamma(1,2,3)


## ----echo=FALSE---------------------------------------------------------------
cat('For example, lambda=1, alpha=2, beta=3: ')
cat('E[X(10)]=',p[1])
cat('Var[X(10)]=',p[2])
cat('lambda*t*E[Y]=',p[3])
cat('lambda*t*E[Y^2]=',p[4])

## -----------------------------------------------------------------------------
x.mean.es<-rep(0,10)
x.var.es<-rep(0,10)
x.mean.th<-rep(0,10)
x.var.th<-rep(0,10)

for (i in 1:10) {
  p<-poissongamma(sample(1:5,1),sample(1:5,1),sample(1:5,1))
  x.mean.es[i]<-p[1]
  x.var.es[i]<-p[2]
  x.mean.th[i]<-p[3]
  x.var.th[i]<-p[4]
}
plot(x.mean.es,x.mean.th,pch=16,xlab=expression(EX(10)),ylab=expression(lambda*t*EY))
abline(1,1)
plot(x.var.es,x.var.th,pch=16,xlab=expression(Var(X(10))),ylab=expression(lambda*t*EY^2))
abline(1,1)
knitr::kable(cbind(x.mean.es,x.mean.th),digits = 1)
knitr::kable(cbind(x.var.es,x.var.th),digits = 1)

## -----------------------------------------------------------------------------
m<-1e4
F<-function(x){#Monte Carlo estimate of the Beta(3, 3)
  u<-runif(m, min = 0, max = x)#u~U(0,x)
  f1<-mean(x*u*u*(1-u)*(1-u))
  v<-runif(m, min = 0, max = 1)#v~U(0,1)
  f2<-mean(v*v*(1-v)*(1-v))
  return(f1/f2)
}
beta_estimated<-rep(0,9)#Monte Carlo estimate of the Beta(3, 3)
beta_true<-rep(0,9)#values returned by the pbeta function
for (i in 1:9) {
  beta_estimated[i]<-F(0.1*i)
  beta_true[i]<-pbeta(0.1*i,3,3)
}
print(round(cbind(seq(0.1,0.9,0.1),beta_estimated,beta_true),5))
plot(beta_estimated,beta_true,pch=16,xlab = "Monte Carlo estimate",ylab = "pbeta function")
abline(0,1,col=2)

## -----------------------------------------------------------------------------
rayleigh<-function(sigma, R=10000, antithetic = TRUE){#inverse transform algorithm 
  u<-runif(R/2)#Uniform(0,1)
  if(!antithetic) 
    v<-runif(R/2) 
  else
    v<-1-u
  x<-(sigma*sqrt(-2*log(1-u))+sigma*sqrt(-2*log(1-v)))/2#(X+X')/2 and (X1+X2)/2
  return(x)
}
sigma<-c(1,3,5,10,20,30)
anti.var<-rep(0,6)#standard deviation using antithetic variables 
notanti.var<-rep(0,6)#standard deviation without using antithetic variables
per<-rep(0,6)#percent reduction in variance
for (i in 1:6) {
  anti.var[i]<-var(rayleigh(sigma = sigma[i]))
  notanti.var[i]<-var(rayleigh(sigma = sigma[i], anti= FALSE))
  per[i]<-100*(notanti.var[i]-anti.var[i])/notanti.var[i]
}
print(round(cbind(sigma,anti.var,notanti.var,per),4))

## -----------------------------------------------------------------------------
m<-100000
theta.hat<-variance<-rep(0,2)
g<-function(x){
  x*x*exp(-x*x/2)*(x>1)/sqrt(2*pi)
}

x<-rnorm(m)#using f1
fg<-g(x)*sqrt(2*pi)/exp(-x*x/2)
theta.hat[1]<-mean(fg)
variance[1]<-var(fg)/m

rayleigh1<-function(sigma, R=10000){#inverse transform algorithm 
  u<-runif(R)#Uniform(0,1)
  x<-sigma*sqrt(-2*log(1-u))#generate random samples from the Rayleigh distribution
  return(x)
}
x<-rayleigh1(sigma = 1, R = m)#using f2
fg<-g(x)/(x*exp(-x*x/2))
theta.hat[2]<-mean(fg)
variance[2]<-var(fg)/m
print(cbind(theta.hat=round(theta.hat,4),variance))

## -----------------------------------------------------------------------------
x<-seq(1.01,4,by=0.01)
plot(x,g(x),type = 'l',lwd=2,xlab = "x",ylab = "f(x)",main = "Importance functions with g(x)",ylim = c(0,0.5))
lines(x,exp(-x*x/2)/sqrt(2*pi),col=2,lwd=2)
lines(x,x*exp(-x*x/2),col=3,lwd=2)
legend("topright",legend = c("g(x)","f1(x)","f2(x)"),lwd = 2,col=c(1,2,3))

plot(x,g(x)*sqrt(2*pi)/exp(-x*x/2),type = 'l',lwd=2,xlab = "x",ylab = "g(x)/f(x)",main = "The ratios g(x)/f(x)",ylim = c(0,6),col=2)
lines(x,g(x)/(x*exp(-x*x/2)),col=3,lwd=2)
legend("topright",legend = c("g(x)/f1(x)","g(x)/f2(x)"),lwd = 2,col=c(2,3))

summary(g(x)*sqrt(2*pi)/exp(-x*x/2))
summary(g(x)/(x*exp(-x*x/2)))

## -----------------------------------------------------------------------------
x<-rayleigh1(sigma = 1, R = m)
fg<-g(x)/(x*exp(-x*x/2))
theta.hat<-mean(fg)
theta.hat

## -----------------------------------------------------------------------------
#卡方分布均值的t-检验
set.seed(1)
n<-20#单次实验的样本量
alpha<-0.05#置信系数
CL<-replicate(1000,expr={
  x<-rchisq(n,df=2)#生成卡方分布随机数
  c(mean(x)-qt(alpha/2,df=n-1)*sd(x)/sqrt(n),mean(x)+qt(alpha/2,df=n-1)*sd(x)/sqrt(n))#计算置信区间
})
s<-0
for (i in 1:1000) {#计算包含真实值的置信区间数量
  if(CL[1,i]>2 & CL[2,i]<2) s<-s+1
}
s/1000#包含真实值的置信区间比例

## -----------------------------------------------------------------------------
set.seed(1)
#正态分布的方差检验
n<-20#单次实验的样本量
alpha<-0.05#置信系数
UCL<-replicate(1000,expr={
  x<-rnorm(n,0,2)#生成正态分布随机样本
  (n-1)*var(x)/qchisq(alpha,df=n-1)#计算置信区间
})
mean(UCL>4)#包含真实值的置信区间比例

#卡方分布的方差检验
n<-20#单次实验的样本量
alpha<-0.05#置信系数
UCL<-replicate(1000,expr={
  x<-rchisq(n,df=2)#生成卡方分布随机样本
  (n-1)*var(x)/qchisq(alpha,df=n-1)#计算置信区间
})
mean(UCL>4)#包含真实值的置信区间比例

## -----------------------------------------------------------------------------
set.seed(1)
n<-20
alpha<-0.05
mu0<-1#真实均值

m<-10000#实验次数
p<-numeric(m)#储存每次实验的p-value
for (i in 1:m) {
  x<-rchisq(n,df=1)#生成卡方分布样本
  ttest<-t.test(x,alternative="two.sided",mu=mu0)#双边t检验
  p[i]<-ttest$p.value
}
p.hat<-mean(p<alpha)#I型错误率
se.hat<-sqrt(p.hat*(1-p.hat)/m)#估计的标准差
print(c(p.hat,se.hat))

## -----------------------------------------------------------------------------
set.seed(1)
n<-20
alpha<-0.05
mu0<-1#真实均值

m<-10000#实验次数
p<-numeric(m)#储存每次实验的p-value
for (i in 1:m) {
  x<-runif(n,min = 0,max = 2)#生成均匀分布样本
  ttest<-t.test(x,alternative="two.sided",mu=mu0)#双边t检验
  p[i]<-ttest$p.value
}
p.hat<-mean(p<alpha)#I型错误率
se.hat<-sqrt(p.hat*(1-p.hat)/m)#估计的标准差
print(c(p.hat,se.hat))

## -----------------------------------------------------------------------------
set.seed(1)
n<-20
alpha<-0.05
mu0<-1#真实均值

m<-10000#实验次数
p<-numeric(m)#储存每次实验的p-value
for (i in 1:m) {
  x<-rexp(n,rate=1)#生成指数分布样本
  ttest<-t.test(x,alternative="two.sided",mu=mu0)#双边t检验
  p[i]<-ttest$p.value
}
p.hat<-mean(p<alpha)#I型错误率
se.hat<-sqrt(p.hat*(1-p.hat)/m)#估计的标准差
print(c(p.hat,se.hat))

## -----------------------------------------------------------------------------
library(knitr)
library(MASS)

## -----------------------------------------------------------------------------
cv.n<-function(d){#样本数量不同时临界值的计算，其中d为数据维度
  n<-c(10,20,30,50,100,500)#样本数量
  cv<-6*qchisq(0.95,df=d*(d+1)*(d+2)/6)/n#不同样本数量下的临界值
  return(cv)
}
d<-3
cv<-cv.n(3)
knitr::kable(rbind(c("n",10,20,30,50,100,500),c("cv",round(cv,4))))

## ----eval=TRUE----------------------------------------------------------------
sk<-function(x){#计算样本X偏度统计量的函数b_{1,d}
  n<-ncol(x)#样本数量
  d<-nrow(x)#样本的数据维度
  xbar<-rowMeans(x)#样本的均值
  xcov<-cov(t(x-xbar))#样本的协方差阵
  x.cen<-x-xbar#样本数据中心化
  return(sum((t(x.cen)%*%solve(xcov)%*%x.cen)^3)/n/n)#输出偏度统计量的函数b_{1,d}
}

## ----eval=TRUE----------------------------------------------------------------
n<-c(10,30,50,100,500)
p.reject<-rep(0,length(n))#储存不同样本数量下拒接原假设的概率
m<-500#每个样本数据下的实验次数
d<-3
for (i in 1:length(n)) {
  sktests<-rep(0,m)
  for (j in 1:m) {
    x<-t(as.matrix(mvrnorm(n=n[i],mu=c(0,0,0),Sigma=matrix(c(3,2,1,2,3,2,1,2,3),3,3))))#生成n[i]个样本
    sktests[j]<-as.integer(abs(sk(x))>=cv[i])
  }
  p.reject[i]<-mean(sktests)#拒绝原假设的概率
}
p.reject
knitr::kable(rbind(c("n",10,30,50,100,500),c("p.reject",round(p.reject,4))))

## ----eval=TRUE----------------------------------------------------------------
alpha<-0.1#置信系数
n<-100#样本数量
m<-500#实验次数
epsilon<-c(seq(0,0.15,0.01),seq(0.15,1,0.05))
N<-length(epsilon)
pwr<-rep(0,N)#功效
cv<-6*qchisq(1-alpha,df=d*(d+1)*(d+2)/6)/n#临界值
s1<-matrix(c(10,2,1,2,10,2,1,2,10),3,3)
s2<-matrix(c(3,2,1,2,3,2,1,2,3),3,3)
x<-matrix(0,nrow=3,ncol=n)

for (j in 1:N) {
  e<-epsilon[j]
  sktests<-rep(0,m)
  for (i in 1:m) {
    for (k in 1:n) {#生成污染的正态分布样本
      t<-sample(c(1,2),1,prob = c(1-e,e))
      if(t==1) x[,k]<-mvrnorm(1,mu=c(0,0,0),Sigma=s1)
      else x[,k]<-mvrnorm(1,mu=c(0,0,0),Sigma=s2)
    }
    sktests[i]<-as.integer(abs(sk(x))>=cv)
  }
  pwr[j]<-mean(sktests)
}

plot(epsilon,pwr,type="b",xlab = bquote(epsilon),ylim = c(0,1),main="Empirical Power")#绘制功效-epsilon图
abline(h=0.1,lty=3)
se<-sqrt(pwr*(1-pwr)/m)#标准差
lines(epsilon,pwr+se,lty=3)
lines(epsilon,pwr-se,lty=3)

## -----------------------------------------------------------------------------
library(bootstrap)
data(scor,package = "bootstrap")
library(boot)
library(knitr)

## -----------------------------------------------------------------------------
set.seed(1)
#样本估计
lambda.hat<-sort(eigen(cov(scor))$values,decreasing = T)#求解协方差矩阵的MLE估计矩阵的特征值并按降序排列
theta.hat<-lambda.hat[1]/sum(lambda.hat)#样本估计

n<-nrow(scor)
B<-2000
theta.b<-rep(0,B)

for (b in 1:B) {#bootstrap
  i<-sample(1:n,size=n,replace = TRUE)
  lambda.b<-sort(eigen(cov(scor[i,]))$values,decreasing = T)
  theta.b[b]<-lambda.b[1]/sum(lambda.b)
}
bias<-mean(theta.b)-theta.hat#偏差
se<-sd(theta.b)#标准差
print(list(theta.hat=theta.hat,bias=bias,se=se))

## -----------------------------------------------------------------------------
#jackknife
theta.jack<-rep(0,n)
for (i in 1:n) {
  lambda.jack<-sort(eigen(cov(scor[-i,]))$values,decreasing = T)
  theta.jack[i]<-lambda.jack[1]/sum(lambda.jack)
}
bias<-(n-1)*(mean(theta.jack)-theta.hat)#偏差
se<-sqrt(n-1)*mean((theta.jack-mean(theta.jack))^2)#标准差
print(list(theta.hat=theta.hat,bias=bias,se=se))

## -----------------------------------------------------------------------------
stat<-function(dat,index){#参数statistic
  lambda<-sort(eigen(cov(dat[index,]))$values,decreasing = T)
  lambda[1]/sum(lambda)
}
boot.obj<-boot(scor,statistic = stat,R=2000)#bootstrap
ci<-boot.ci(boot.obj,type=c("perc","bca"))#计算percentile和BCa置信区间
knitr::kable(rbind(c("precent",ci$perc[4],ci$perc[5]),c("BCa",ci$bca[4],ci$bca[5])),digits = 4)

## ----eval=TRUE----------------------------------------------------------------
sk<-function(x,index){
  xbar<-mean(x[index])
  m3<-mean((x[index]-xbar)^3)
  m2<-mean((x[index]-xbar)^2)
  m3/m2^1.5
}

## ----eval=TRUE----------------------------------------------------------------
set.seed(1)
m<-1000#蒙特卡洛实验次数
x<-rnorm(10000)#正态样本
sk.x<-mean((x-mean(x))^3)/(mean((x-mean(x))^2))^1.5#计算样本的偏度系数
sk.x
sktests<-matrix(0,nrow=6,ncol=m)
sk.true<-matrix(0,nrow=6,ncol=m)
for (j in 1:m) {#进行蒙特卡洛实验
  x<-rnorm(1000)
  #bootstrap
  boot.obj<-boot(x,statistic = sk,R=100)
  b<-boot.ci(boot.obj,type = c("norm", "basic", "perc"))
  sktests[1,j]<-as.integer(sk.x<=b$normal[2])#样本偏度系数normal检验落入左侧
  sktests[2,j]<-as.integer(sk.x>=b$normal[3])#样本偏度系数normal检验落入右侧
  sktests[3,j]<-as.integer(sk.x<=b$basic[4])#样本偏度系数basic检验落入左侧
  sktests[4,j]<-as.integer(sk.x>=b$basic[5])#样本偏度系数basic检验落入右侧
  sktests[5,j]<-as.integer(sk.x<=b$percent[4])#样本偏度系数percent检验落入左侧
  sktests[6,j]<-as.integer(sk.x>=b$percent[5])#样本偏度系数percent检验落入右侧
  
  sk.true[1,j]<-as.integer(0<=b$normal[2])#真实偏度系数normal检验落入左侧
  sk.true[2,j]<-as.integer(0>=b$normal[3])#真实偏度系数normal检验落入右侧
  sk.true[3,j]<-as.integer(0<=b$basic[4])#真实偏度系数basic检验落入左侧
  sk.true[4,j]<-as.integer(0>=b$basic[5])#真实偏度系数basic检验落入右侧
  sk.true[5,j]<-as.integer(0<=b$percent[4])#真实偏度系数percent检验落入左侧
  sk.true[6,j]<-as.integer(0>=b$percent[5])#真实偏度系数percent检验落入右侧
}
samplesk<-rowMeans(sktests)
truesk<-rowMeans(sk.true)
cat("样本偏度系数")
knitr::kable(rbind(c("type","miss on the left","miss on the right","coverage rate"),c("normal",samplesk[1:2],1-sum(samplesk[1:2])),c("basic",samplesk[3:4],1-sum(samplesk[3:4])),c("percet",samplesk[5:6],1-sum(samplesk[5:6]))))
cat("真实偏度系数")
knitr::kable(rbind(c("type","miss on the left","miss on the right","coverage rate"),c("normal",truesk[1:2],1-sum(truesk[1:2])),c("basic",truesk[3:4],1-sum(truesk[3:4])),c("percet",truesk[5:6],1-sum(truesk[5:6]))))

## ----eval=TRUE----------------------------------------------------------------
set.seed(1)
m<-1000#蒙特卡洛实验次数
x<-rchisq(10000,df=5)#正态样本
sk.x<-mean((x-mean(x))^3)/(mean((x-mean(x))^2))^1.5#计算样本的偏度系数
sk.x
sktests<-matrix(0,nrow=6,ncol=m)
sk.true<-matrix(0,nrow=6,ncol=m)
for (j in 1:m) {#进行蒙特卡洛实验
  x<-rchisq(1000,df=5)
  #bootstrap
  boot.obj<-boot(x,statistic = sk,R=100)
  b<-boot.ci(boot.obj,type = c("norm", "basic", "perc"))
  sktests[1,j]<-as.integer(sk.x<=b$normal[2])#样本偏度系数normal检验落入左侧
  sktests[2,j]<-as.integer(sk.x>=b$normal[3])#样本偏度系数normal检验落入右侧
  sktests[3,j]<-as.integer(sk.x<=b$basic[4])#样本偏度系数basic检验落入左侧
  sktests[4,j]<-as.integer(sk.x>=b$basic[5])#样本偏度系数basic检验落入右侧
  sktests[5,j]<-as.integer(sk.x<=b$percent[4])#样本偏度系数percent检验落入左侧
  sktests[6,j]<-as.integer(sk.x>=b$percent[5])#样本偏度系数percent检验落入右侧
  
  sk.true[1,j]<-as.integer(1.265<=b$normal[2])#真实偏度系数normal检验落入左侧
  sk.true[2,j]<-as.integer(1.265>=b$normal[3])#真实偏度系数normal检验落入右侧
  sk.true[3,j]<-as.integer(1.265<=b$basic[4])#真实偏度系数basic检验落入左侧
  sk.true[4,j]<-as.integer(1.265>=b$basic[5])#真实偏度系数basic检验落入右侧
  sk.true[5,j]<-as.integer(1.265<=b$percent[4])#真实偏度系数percent检验落入左侧
  sk.true[6,j]<-as.integer(1.265>=b$percent[5])#真实偏度系数percent检验落入右侧
}
samplesk<-rowMeans(sktests)
truesk<-rowMeans(sk.true)
cat("样本偏度系数")
knitr::kable(rbind(c("type","miss on the left","miss on the right","coverage rate"),c("normal",samplesk[1:2],1-sum(samplesk[1:2])),c("basic",samplesk[3:4],1-sum(samplesk[3:4])),c("percet",samplesk[5:6],1-sum(samplesk[5:6]))))
cat("真实偏度系数")
knitr::kable(rbind(c("type","miss on the left","miss on the right","coverage rate"),c("normal",truesk[1:2],1-sum(truesk[1:2])),c("basic",truesk[3:4],1-sum(truesk[3:4])),c("percet",truesk[5:6],1-sum(truesk[5:6]))))

## -----------------------------------------------------------------------------
library(knitr)
library(MASS)
library(boot)
library(RANN)
library(energy)
library(Ball)

## -----------------------------------------------------------------------------
set.seed(1)
x<-mvrnorm(100,c(0,0),Sigma = matrix(c(2,1,1,2),2,2))
y<-mvrnorm(100,c(0,0),Sigma = diag(2))

## -----------------------------------------------------------------------------
spear<-function(z,ix){
  #dims是x和y的维度
  x<-z[,1:2]#x不改变
  y<-z[ix,3:4]#重排y的行
  return(cor.test(x,y,method = "spearman")$p.value)
}

z<-cbind(x,y)
boot.obj<-boot(data=z,statistic = spear,R=999,sim="permutation")
tb<-c(boot.obj$t0,boot.obj$t)
mean(tb<=boot.obj$t0)

## -----------------------------------------------------------------------------
x<-mvrnorm(100,c(0,0),Sigma = diag(2))
y<-x

## -----------------------------------------------------------------------------
spear<-function(z,ix){
  #dims是x和y的维度
  x<-z[,1:2]#x不改变
  y<-z[ix,3:4]#重排y的行
  return(cor.test(x,y,method = "spearman")$p.value)
}
z<-cbind(x,y)
boot.obj<-boot(data=z,statistic = spear,R=999,sim="permutation")
tb<-c(boot.obj$t0,boot.obj$t)
mean(tb<=boot.obj$t0)

## ----eval=TRUE----------------------------------------------------------------
Tn<-function(z,ix,sizes,k){#进行NN的统计量计算
  n1<-sizes[1]
  n2<-sizes[2]
  n<-n1+n2
  if(is.vector(z)) z<-data.frame(z,0)#确保z不是向量以保证以下代码的正确性
  z<-z[ix,]
  NN<-nn2(data=z,k=k+1)#NN
  block1<-NN$nn.idx[1:n1,-1]#寻找前n1个邻居
  block2<-NN$nn.idx[(n1+1):n,-1]#寻找后n2个邻居
  i1<-sum(block1<n1+0.5)
  i2<-sum(block2>n1+0.5)
  (i1+i2)/(k*n)
}

eqdist.nn<-function(z,sizes,k){
  boot.obj<-boot(data=z,statistic = Tn,R=R,sim="permutation",sizes=sizes,k=k)#重排法
  ts<-c(boot.obj$t0,boot.obj$t)
  p.value<-mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

## ----eval=TRUE----------------------------------------------------------------
m<-200
k<-3
p<-2#数据维度
set.seed(1)
n1<-n2<-50#样本数量
R<-100
n<-n1+n2
N<-c(n1,n2)

p.values<-matrix(0,m,3)
for (i in 1:m) {
  x<-mvrnorm(n1,c(0,0),diag(2))
  y<-mvrnorm(n2,c(0,0),diag(2,2))
  z<-rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#ENERGY
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*1)$p.value#BALL
}
alpha<-0.1
pow<-colMeans(p.values<alpha)
pow
alpha<-0.05
pow<-colMeans(p.values<alpha)
pow

## ----eval=TRUE----------------------------------------------------------------
m<-200
k<-3
p<-2#数据维度
set.seed(1)
n1<-n2<-50#样本数量
R<-100
n<-n1+n2
N<-c(n1,n2)

p.values<-matrix(0,m,3)
for (i in 1:m) {
  x<-mvrnorm(n1,c(0,0),diag(2))
  y<-mvrnorm(n2,c(1,1),diag(2,2))
  z<-rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#ENERGY
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*1)$p.value#BALL
}

alpha<-0.05
pow<-colMeans(p.values<alpha)
pow
alpha<-0.005
pow<-colMeans(p.values<alpha)
pow

## ----eval=TRUE----------------------------------------------------------------
m<-200
k<-3
p<-2#数据维度
set.seed(1)
n1<-n2<-50#样本数量
R<-100
n<-n1+n2
N<-c(n1,n2)

p.values<-matrix(0,m,3)
for (i in 1:m) {
  x<-as.matrix(rt(n1,df=1))
  s1<-rbinom(1,n2,0.3)
  s2<-n2-s1
  y1<-rnorm(s1,0,1)
  y2<-rnorm(s2,3,2)
  y<-c(y1,y2)
  y<-matrix(y)
  z<-rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#ENERGY
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*1)$p.value#BALL
}
alpha<-0.05
pow<-colMeans(p.values<alpha)
pow
alpha<-0.01
pow<-colMeans(p.values<alpha)
pow

## ----eval=TRUE----------------------------------------------------------------
m<-100
k<-3
p<-2#数据维度
set.seed(1)
n1<-10
n2<-90#样本数量
R<-100
n<-n1+n2
N<-c(n1,n2)

p.values<-matrix(0,m,3)
for (i in 1:m) {
  x<-mvrnorm(n1,c(0,0),diag(2))
  y<-mvrnorm(n2,c(0,0),diag(2,2))
  z<-rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#ENERGY
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*1)$p.value#BALL
}
alpha<-0.1
pow<-colMeans(p.values<alpha)
pow
alpha<-0.05
pow<-colMeans(p.values<alpha)
pow


## ----eval=TRUE----------------------------------------------------------------
m<-200
k<-3
p<-2#数据维度
set.seed(1)
n1<-10
n2<-90#样本数量
R<-100
n<-n1+n2
N<-c(n1,n2)

p.values<-matrix(0,m,3)
for (i in 1:m) {
  x<-as.matrix(rt(n1,df=1))
  s1<-rbinom(1,n2,0.3)
  s2<-n2-s1
  y1<-rnorm(s1,0,1)
  y2<-rnorm(s2,3,2)
  y<-c(y1,y2)
  y<-matrix(y)
  z<-rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value#NN
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value#ENERGY
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,seed=i*1)$p.value#BALL
}
alpha<-0.1
pow<-colMeans(p.values<alpha)
pow
alpha<-0.05
pow<-colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
library(knitr)
library(coda)

## -----------------------------------------------------------------------------
set.seed(123)
#标准柯西分布的密度函数
f<-function(x){
  return(1/(pi*(1+x^2)))
}

#生成标准柯西分布的马尔科夫链
cauchy.chain<-function(sigma,N,X1){
  #sigma为提议分布的标准差
  #N为链的长度
  #初始值为X1
  x<-rep(0,N)
  x[1]<-X1
  u<-runif(N)
  
  for (i in 2:N) {
    xt<-x[i-1]
    y<-rnorm(1,xt,sigma)
    num<-f(y)*dnorm(xt,y,sigma)
    den<-f(xt)*dnorm(y,xt,sigma)
    r<-num/den
    if(u[i]<=r) x[i]<-y
    else x[i]<-xt
  }
  return(x)
}

#Gelman-Rubin method
gelman.rubin<-function(psi){
  psi<-as.matrix(psi)
  n<-ncol(psi)
  k<-nrow(psi)
  
  psi.means<-rowMeans(psi)
  B<-n*var(psi.means)
  psi.w<-apply(psi,1,"var")
  W<-mean(psi.w)
  v.hat<-W*(n-1)/n+(B/n)
  r.hat<-v.hat/W
  return(r.hat)
}

sigma<-1#提议分布的标准差
n<-11000#链的长度
b<-1000#预烧期值

X<-cauchy.chain(sigma,n,0)
X<-X[-(1:b)]

## ----eval=TRUE----------------------------------------------------------------
#实验生成的随机数与真实Cauchy分布分位数的比较
a<-0.1*c(1:9)
knitr::kable(rbind(c("generated observations",round(quantile(X,a),4)),c("standard Cauchy distribution",round(qcauchy(a),4))))
plot(quantile(X,a),qcauchy(a),type="p",xlab = "generated observations",ylab = "standard Cauchy distribution")
abline(a=0,b=1)

## ----eval=TRUE----------------------------------------------------------------
set.seed(1)
sigma<-1#提议分布的标准差
k<-4#生成的马尔可夫链数目
n<-20000#链的长度
b<-1000#预烧期值

#生成链
X<-matrix(0,nrow=k,ncol=n-b)
for (i in 1:k) {
  X[i,]<-cauchy.chain(sigma,n,0)[-c(1:b)]#去掉前1000个迭代值
}
psi<-t(apply(X,1,cumsum))
for (i in 1:nrow(psi)) {
  psi[i,]<-psi[i,]/(1:ncol(psi))
}

rhat<-rep(0,n-b)
for(j in 5000:(n-b)){
  rhat[j]<-gelman.rubin(psi[,1:j])#从链长度为5000时开始计算
}
rhat<-rhat[5000:(n-b)]
plot(5000:(n-b),rhat,type="l")
abline(h=1.2,lty=2)
which(rhat<=1.2)[1]+5000

## -----------------------------------------------------------------------------
set.seed(1)
sigma<-3#提议分布的标准差
k<-4#生成的马尔可夫链数目
n<-20000#链的长度
b<-1000#预烧期值

#生成链
X<-matrix(0,nrow=k,ncol=n-b)
for (i in 1:k) {
  X[i,]<-cauchy.chain(sigma,n,0)[-c(1:b)]#去掉前1000个迭代值
}
psi<-t(apply(X,1,cumsum))
for (i in 1:nrow(psi)) {
  psi[i,]<-psi[i,]/(1:ncol(psi))
}

rhat<-rep(0,n-b)
for(j in 5000:(n-b)){
  rhat[j]<-gelman.rubin(psi[,1:j])#从链长度为5000时开始计算
}
rhat<-rhat[5000:(n-b)]
plot(5000:(n-b),rhat,type="l")
abline(h=1.2,lty=2)
which(rhat<=1.2)[1]+5000

## -----------------------------------------------------------------------------
N<-10000#链的长度
burn<-1000#预烧期值


#生成二元的马尔可夫链
chain<-function(n,a,b,N,X1){
  set.seed(1)
  #sigma为提议分布的标准差
  #N为链的长度
  #初始值为X1
  X<-matrix(0,N,2)#马尔可夫链
  X[1,]<-X1
  
  
  for (i in 2:N) {
    y<-X[i-1,2]
    x<-X[i,1]<-rbinom(1,n,prob = y)
    X[i,2]<-rbeta(1,shape1=x+a,shape2=n-x+b)
  }
  return(X[(burn+1):N,])
}

## -----------------------------------------------------------------------------
set.seed(1)
n<-19
a<-1
b<-1
#更改初始值
x1<-sample(1:n,1)
y1<-runif(1)
X1<-c(x1,y1)
x2<-sample(1:n,1)
y2<-runif(1)
X2<-c(x2,y2)
x3<-sample(1:n,1)
y3<-runif(1)
X3<-c(x3,y3)

#生成三个初始值不同的马尔可夫链
X<-list(chain(n,a,b,N,X1),chain(n,a,b,N,X2),chain(n,a,b,N,X3))
#Gelman-Rubin诊断
mcmc<-mcmc.list(mcmc(X[[1]]),mcmc(X[[2]]),mcmc(X[[3]]))
gelman.diag(mcmc,transform = T)
gelman.plot(mcmc,transform = T)

## -----------------------------------------------------------------------------
set.seed(1)
n<-5
a<-1
b<-1
#更改初始值
x1<-sample(1:n,1)
y1<-runif(1)
X1<-c(x1,y1)
x2<-sample(1:n,1)
y2<-runif(1)
X2<-c(x2,y2)
x3<-sample(1:n,1)
y3<-runif(1)
X3<-c(x3,y3)

#生成三个初始值不同的马尔可夫链
X<-list(chain(n,a,b,N,X1),chain(n,a,b,N,X2),chain(n,a,b,N,X3))
#Gelman-Rubin诊断
mcmc<-mcmc.list(mcmc(X[[1]]),mcmc(X[[2]]),mcmc(X[[3]]))
gelman.diag(mcmc,transform = T)
gelman.plot(mcmc,transform = T)

## -----------------------------------------------------------------------------
#计算第k项的值
k.term<-function(k,a,d){
  #k>=0,d>0整数，a是一个d维向量
  (-1)^k*exp(-log(factorial(k)) - k*log(2) + (k+1)*log(sum(a^2)) - log(2*k+1) - log(2*k+2) + lgamma((d+1)/2) + lgamma(k+1.5) - lgamma(k+(d/2)+1))
}

## -----------------------------------------------------------------------------
#计算求和
k.sum<-function(a,d){
  #d>0整数，a是一个d维向量
  s<-0
  kterm<-1
  k<-0
  while (k<=1000 & abs(kterm)>10e-10) {#当迭代次数超过1000或者单项绝对值小于10e-10停止
    kterm<-k.term(k,a,d)
    s<-s+kterm
    k<-k+1
  }
  return(s)
}

## -----------------------------------------------------------------------------
a<-c(1,2)
d<-2
k.sum(a,d)

## -----------------------------------------------------------------------------
#积分函数
int<-function(a,k){
  ck<-sqrt(a*a*k/(k+1-a*a))
  integrate(function(u) {(1+u*u/k)^(-(k+1)/2)},lower = 0,upper = ck)$value
}

Sk<-function(a,k){#计算
  2*exp(lgamma((k+1)/2)-log(sqrt(pi*k))-lgamma(k/2))*int(a,k)
}

k<-c(4:25,100,500,1000)
a<-seq(0,3,0.1)
y<-matrix(0,nrow=4,ncol = length(a))
for (i in 1:length(a)) {
  y[1,i]<-Sk(a[i],k[22])-Sk(a[i],k[22]-1)
  y[2,i]<-Sk(a[i],k[23])-Sk(a[i],k[23]-1)
  y[3,i]<-Sk(a[i],k[24])-Sk(a[i],k[24]-1)
  y[4,i]<-Sk(a[i],k[25])-Sk(a[i],k[25]-1)
}
#通过图来判断根所在的区间

plot(a,y[1,],type = "l",xlab = "a",ylab = "S_x-S_{k-1}",main = "k=25")
abline(h=0,col=2)

plot(a,y[2,],type = "l",xlab = "a",ylab = "S_x-S_{k-1}",main = "k=100")
abline(h=0,col=2)

plot(a,y[3,],type = "l",xlab = "a",ylab = "S_x-S_{k-1}",main = "k=500")
abline(h=0,col=2)

plot(a,y[4,],type = "l",xlab = "a",ylab = "S_x-S_{k-1}",main = "k=1000")
abline(h=0,col=2)

## -----------------------------------------------------------------------------
Ak<-function(k){
  uniroot(function(a) {Sk(a,k)-Sk(a,k-1)},interval=c(1,2))$root
}
root1<-rep(0,length(k))
for (i in 1:length(k)) {
  root1[i]<-Ak(k[i])
}
knitr::kable(rbind(c("k",k),c("root",round(root1,3))))

## -----------------------------------------------------------------------------
Sk<-function(a,k){
  critical<-sqrt(a*a*k/(k+1-a*a))
  1-pt(critical,df=k)
}
k<-c(4:25,100,500,1000)
a<-seq(0,3,0.1)
#通过图来判断根所在的区间

plot(a,Sk(a,k[22])-Sk(seq(0,3,0.1),k[22]-1),type = "l",xlab = "a",ylab = "S_x-S_{k-1}",main = "k=25")
abline(h=0,col=2)

plot(a,Sk(a,k[23])-Sk(seq(0,3,0.1),k[23]-1),type = "l",xlab = "a",ylab = "S_x-S_{k-1}",main = "k=100")
abline(h=0,col=2)

plot(a,Sk(a,k[24])-Sk(seq(0,3,0.1),k[24]-1),type = "l",xlab = "a",ylab = "S_x-S_{k-1}",main = "k=500")
abline(h=0,col=2)

plot(a,Sk(a,k[25])-Sk(seq(0,3,0.1),k[25]-1),type = "l",xlab = "a",ylab = "S_x-S_{k-1}",main = "k=1000")
abline(h=0,col=2)


## -----------------------------------------------------------------------------
Ak<-function(k){#求解A(k)
  uniroot(function(a) {Sk(a,k)-Sk(a,k-1)},interval=c(1,2))$root
}
root2<-rep(0,length(k))
for (i in 1:length(k)) {
  root2[i]<-Ak(k[i])
}
knitr::kable(rbind(c("k",k),c("root",round(root2,3))))

## -----------------------------------------------------------------------------
knitr::kable(rbind(c("k",k),c("11.4",round(root2,3)),c("11.5",round(root1,3)),c("差",round(root1-root2,3))))

## -----------------------------------------------------------------------------
y<-c(0.54,0.48,0.33,0.43,0.91,0.21,0.85,1.00,1.00,1.00)
l<-function(lambda){#似然函数
  -7*log(lambda)-sum(y[1:7])/lambda-3/lambda
}
optimize(l,interval = c(0.01,3),maximum = T)$maximum

## -----------------------------------------------------------------------------
#E-M算法
EM<-function(y,max.it=10000,eps=1e-5){
  i<-1
  lambda1<-1
  lambda2<-2
  while (abs(lambda1-lambda2)>eps & i<max.it) {
    lambda1<-lambda2
    Q<-function(lambda){#E步
      -7*log(lambda)-sum(y[1:7])/lambda-3*(1+lambda1)/lambda
    }
    lambda2<-optimize(Q,interval = c(0.01,3),maximum = T)$maximum#M步
    i<-i+1
  }
  return(lambda2)
}
y<-c(0.54,0.48,0.33,0.43,0.91,0.21,0.85,1.00,1.00,1.00)
EM(y)

## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
data(mtcars)
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
rsq <- function(mod) summary(mod)$r.squared
lapply(formulas,function(x){
  mod<-lm(x,data = mtcars)
  rsq(mod)
})

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
rsq <- function(mod) summary(mod)$r.squared
lapply(bootstraps,function(x){
  mod<-lm(mpg ~ disp,data=x)
  rsq(mod)
})

## -----------------------------------------------------------------------------
x1<-1:20
x2<-rnorm(20)
x3<-rep(5,20)
x<-data.frame(x1,x2,x3)#生成一个全是数据的dataframe
vapply(x,sd,numeric(1))#求解dataframe每一列数据的标准差
print(c(sd(x1),sd(x2),sd(x3)))#验证结果正确

## -----------------------------------------------------------------------------
x1<-rnorm(5)
x2<-c("a","b","c","d","e")
x3<-rep(5,5)
x4<-c("1","2","3","4","5")
x<-data.frame(x1,x2,x3,x4)#生成一个混合的dataframe
y<-vapply(x,is.numeric,logical(1))#判断全是数据的列
vapply(x[y],sd,numeric(1))#求数据列的标准差

## -----------------------------------------------------------------------------
library(parallel)
cores<-detectCores()
cores
cluster<-makePSOCKcluster(cores)

boot_rsq <- function(i) {
  summary(lm(mpg ~ wt + disp, data = mtcars[sample(nrow(mtcars), rep = T), ]))$r.square
}

mcsapply<-function(x,f,...){
  res<-parLapply(cluster,x,f,...)
  simplify2array(res)
}

#简单检测正确性
sapply(1:10,sqrt)
mcsapply(1:10,sqrt)
parSapply(cluster,1:10,sqrt)
#比较运行时间
system.time(sapply(1:500,boot_rsq))
system.time(mcsapply(1:500,boot_rsq))
system.time(parSapply(cluster,1:500,boot_rsq))

## -----------------------------------------------------------------------------
vapply2<-function(x,f,f.value,...){
  
  out<-matrix(rep(f.value,length(x)),nrow=length(x))
  for (i in sample(seq_along(x))) {#打乱计算顺序
    res<-f(unlist(x[i]),...)
    stopifnot(
      length(res)==length(f.value),
      typeof(res)==typeof(f.value)
    )
    out[i,]<-res

  }
  out
}

x1<-1:20
x2<-rnorm(20)
x3<-rep(5,20)
x<-data.frame(x1,x2,x3)#生成一个全是数据的dataframe
vapply(x,sd,numeric(1))#求解dataframe每一列数据的标准差
vapply2(x,sd,numeric(1))

x1<-rnorm(5)
x2<-c("a","b","c","d","e")
x3<-rep(5,5)
x4<-c("1","2","3","4","5")
x<-data.frame(x1,x2,x3,x4)#生成一个混合的dataframe
vapply(x,is.numeric,logical(1))#判断全是数据的列
vapply2(x,is.numeric,logical(1))

## -----------------------------------------------------------------------------
cores<-detectCores()
cores
cluster<-makePSOCKcluster(cores)

mcvapply<-function(x,f,f.value,...){
  res<-mcsapply(x,f,...)
  for (i in sample(seq_along(x)))
    stopifnot(
      length(res[i])==length(f.value),
      typeof(res[i])==typeof(f.value)
    )
  res
}

x1<-1:20
x2<-rnorm(20)
x3<-rep(5,20)
x<-data.frame(x1,x2,x3)#生成一个全是数据的dataframe
mcvapply(x,sd,numeric(1))#求解dataframe每一列数据的标准差


x1<-rnorm(5)
x2<-c("a","b","c","d","e")
x3<-rep(5,5)
x4<-c("1","2","3","4","5")
x<-data.frame(x1,x2,x3,x4)#生成一个混合的dataframe
mcvapply(x,is.numeric,logical(1))#判断全是数据的列

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

## -----------------------------------------------------------------------------
set.seed(123)

##Gibbs抽样的Cpp函数
cppFunction('NumericMatrix gibbsC(int n, double a, double b, int N, double x1, double y1) {
  NumericMatrix X(N,2);
  X(0,0)=x1;
  X(0,1)=y1;
  int i;
  double x;
  double y;
  
  for(i=1;i<N;i++){
    y=X(i-1,1);
    X(i,0)=rbinom(1,n,y)[0];
    x=X(i,0);
    X(i,1)=rbeta(1,x+a,n-x+b)[0];
  }
  return(X);
}')

##Gibbs抽样的R函数
gibbsR<-function(n,a,b,N,x1,y1){
  set.seed(1)
  #sigma为提议分布的标准差
  #N为链的长度
  #初始值为X1
  X<-matrix(0,N,2)#马尔可夫链
  X[1,]<-c(x1,y1)
  
  for (i in 2:N) {
    y<-X[i-1,2]
    x<-X[i,1]<-rbinom(1,n,prob = y)
    X[i,2]<-rbeta(1,shape1=x+a,shape2=n-x+b)
  }
  return(X)
}

## -----------------------------------------------------------------------------
#通过qq图比较两种语言写的函数生成的随机数
set.seed(123)
X<-gibbsC(10,2,3,1000,2,0.5)[-(1:500),]
Y<-gibbsR(10,2,3,1000,2,0.5)[-(1:500),]
qqplot(X[,1],Y[,1],xlab = "x of Cpp",ylab = "x of R")
abline(a=0,b=1)
qqplot(X[,2],Y[,2],xlab = "y of Cpp",ylab = "y of R")
abline(a=0,b=1)

## -----------------------------------------------------------------------------
set.seed(123)
microbenchmark(gibbsC=gibbsC(10,2,3,1000,2,0.5),gibbsR=gibbsR(10,2,3,1000,2,0.5))

## -----------------------------------------------------------------------------
#通过qq图比较两种语言写的函数生成的随机数
set.seed(123)
X<-gibbsC(20,7,3,10000,10,0.1)[-(1:500),]
Y<-gibbsR(20,7,3,10000,10,0.1)[-(1:500),]
qqplot(X[,1],Y[,1],xlab = "x of Cpp",ylab = "x of R")
abline(a=0,b=1)
qqplot(X[,2],Y[,2],xlab = "y of Cpp",ylab = "y of R")
abline(a=0,b=1)

## -----------------------------------------------------------------------------
set.seed(123)
microbenchmark(gibbsC=gibbsC(20,7,3,10000,10,0.1),gibbsR=gibbsR(20,7,3,10000,10,0.1))


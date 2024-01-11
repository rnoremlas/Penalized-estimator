rm(list=ls())

library(multiColl)
source("03_functions.txt")
source("04_functions.txt")

############################################

datos = read.table("01_datos_multicol_alta.txt", header=T, sep=";")

Y = datos[,1]
X = as.matrix(datos[,-1])

CVs(X)
VIF(X)
CNs(X)
RdetR(X)

####################################################################

reg = lm(Y~X+0)
summary(reg)

##
 
h = 1
k = 0.01
betahk = BetaKH(Y, X, h, k)
betahk
varBetaKH(Y, X, h, k)

##

p = ncol(X)
alfa = matrix(mean(Y), p, 1)
for (i in 2:p) alfa[i] = cov(Y, X[,i])/var(X[,i])
BetaKH_bootstrap(Y, X, h, k, alfa)

##

BA = 1 - sum((Y-X%*%betahk[[2]])^2)/crossprod(Y)
BA

####################################################################
# bootstrap

ite = 10000

n = nrow(X)
p = ncol(X)

boots = matrix(, ite, p)
BAs = matrix(, ite, 1)

for (i in 1:ite){
  s = sample(1:n, n, replace=T)
  l = 1
  Xb = matrix(, n, p)
  Yb = matrix(, n, 1)
  for (j in s)
  {
    Xb[l,] = X[j,]
    Yb[l] = Y[j]
    l = l + 1
  }
  boots[i,] = as.matrix(BetaKH_bootstrap(Yb, Xb, h, k, alfa)[[2]])
  BAs[i] = 1 - sum((Yb-Xb%*%boots[i,])^2)/crossprod(Yb)
}

########## BETA

# punctual estimation

  beta = colMeans(boots)
  beta
  for (i in 1:p){
    print(sd(boots[,i]))
  }

# interval estimation
  
  # type 1

  int1low = matrix(, p, 1)
  int1top = matrix(, p, 1)
  for (i in 1:p){
    int1low[i] = quantile(boots[,i], prob=0.025)
    int1top[i] = quantile(boots[,i], prob=0.975)
  }
  int1 = cbind(int1low, int1top)
  int1

  # type 2

  int2low = array(,p)
  int2top = array(,p)
  for (i in 1:p){
    int2low[i] = mean(boots[,i])-1.96*sd(boots[,i])
    int2top[i] = mean(boots[,i])+1.96*sd(boots[,i])
  }    
  int2 = cbind(int2low, int2top)
  int2

########## BA

  mean(BAs)
  sd(BAs)
  intBA1 = c(quantile(BAs, prob=0.025), quantile(BAs, prob=0.975))
  intBA1
  intBA2 = c(mean(BAs)-1.96*sd(BAs), mean(BAs)+1.96*sd(BAs))
  intBA2

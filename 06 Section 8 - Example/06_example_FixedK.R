rm(list=ls()) 

source("02_functions.txt")
source("03_functions.txt")
source("04_functions.txt")
source("05_functions.txt")
source("06_functions.txt")
library(multiColl)

datos = read.table("06_Wissel.txt", header=T, sep=";")
head(datos)

attach(datos)

  cte = rep(1,length(D))
  X = cbind(cte, C, I , CP)

  CVs(X)
  VIF(X)
  CNs(X)
  RdetR(X)

  n = nrow(X)
  p = ncol(X)

  ##

  reg = lm(D~C+I+CP)
  summary(reg)
  beta = matrix(reg$coef, p, 1)
  beta
  sigma = as.numeric(summary(reg)[6])
  sigma
  
  ##

  alfa = matrix(mean(D), p, 1)
  alfa[2] = cov(D,C)/var(C)
  alfa[3] = cov(D,I)/var(I)
  alfa[4] = cov(D,CP)/var(CP)
  alfa 
  
  #####################################################
  # for a fixed k
  
  h = 1
  k = 100
  
  PenalizedEstimation(D, X, h, k, alfa, sigma, beta, ite = 10000, tol = 0.01, setseed = 1)
  
detach(datos)







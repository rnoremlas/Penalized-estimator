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
  
  ##
  
  reg.simple2 = lm(D~C)
  summary(reg.simple2) # match with alfa[2]
  reg.simple3 = lm(D~I)
  summary(reg.simple3)# match with alfa[3]
  reg.simple4 = lm(D~CP)
  summary(reg.simple4) # match with alfa[4]

  #####################################################
  # traces

  start = 0
  leap = 0.01
  stop = 1
  
  h=0 # cresta
  
  ##
  
  output = TracesHK(D, X, alfa, h, sigma, beta, start, leap, stop, graph=1)
  # output:
    kas = output[,1] # position 1: k
    # position 2 to p+1: betashk
    # position p+2: norm
    norms.rate = output[, p+3] # position p+3: rate norm
    MSEs = output[, p+4] # position p+4: MSE
    NCs = output[, p+5] # position p+5: NC
    VIFs = output[, (p+6):(2*p+4)] # position p+6 to 2p+4: VIFs
    detRs = output[, 2*p+5] # position 2p+5: detRs
ite=100
    perturb_ite(D, X, h, discretization=0.02, alfa, ite, mu=5, dv=5, tol=0.01)
    
  #####################################################
  # choice of k
  
  firstV(norms.rate, threshold=0.1, kas)
  
  minMSE(MSEs, kas)
  
  firstV(NCs, threshold=30, kas)
  firstV(NCs, threshold=20, kas)
  firstV(NCs, threshold=10, kas)
  
  firstV(VIFs[,1], threshold=10, kas)
  firstV(VIFs[,2], threshold=10, kas)
  firstV(VIFs[,3], threshold=10, kas)
  firstVIF(VIFs, threshold=10, leap, stop)
  
  firstV(detRs, threshold=0.1, kas, minor = F)
  
  #####################################################
  # for a fixed k
  
  ite = 10000
  setseed = 1
  
  fixed_kas = c(0, 0.01, 0.02, 0.04, 0.08) # k = 0 OLS, k = 0.01 NC < 20, k = 0.02 MSE min, k = 0.04 NC < 10, k = 0.08 max VIF < 10
  table = matrix(, 26, 2+2*length(fixed_kas))
  i = 1
  for (k in fixed_kas){
    PenalizedEstimation(D, X, h, k, alfa, sigma, beta, ite, tol = 0.01, setseed = 1)
    output = read.table(paste("06_output_h=", h, "_k=", k, ".txt", sep=""), sep="&")
    table[, 2*i] = round(output[,2], digits=4)
    table[, 2*i+1] = round(output[,3], digits=4)
    i = i + 1
  }
  table[,1] = output[,1]
  table[,2+2*length(fixed_kas)] = output[,4]
  #table
  write.table(table, paste("06_Table_LaTeX_h=", h, ".txt", sep=""), sep=" & ", col.names = F, row.names = F) # open and delete "" for LaTeX document
  
detach(datos)

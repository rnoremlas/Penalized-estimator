rm(list=ls()) 

source("02_functions.txt")
source("03_functions.txt")

###########################################################################

datos = read.table("01_datos_multicol_alta.txt", header=T, sep=";")

Y = datos[,1]
X = as.matrix(datos[,-1])

n = nrow(X)
p = ncol(X)

##

library(multiColl)
  CVs(X)
  RdetR(X)
  VIF(X)
  CNs(X)

##

reg = lm(Y~X+0)
beta = matrix(reg$coefficients, 3, 1)
sigma = summary(reg)[[6]]
alfa = matrix(mean(Y), p, 1)
for (q in 2:p) alfa[q] = cov(Y, X[,q])/var(X[,q])  

###########################################################################
# traces

start = 0
leap = 0.01
stop = 100

##

discretizacion = seq(start, stop, leap)
l = length(discretizacion)

norms = numeric()
MSEs = array(,length(discretizacion))
NCs = array(,length(discretizacion))
VIFs = matrix(,length(discretizacion), p-1)
detRs = array(,length(discretizacion))
i = 1
for (k in discretizacion) {
  betaskh = BetaKH(Y, X, h=1, k)[[2]]
  norms[i] = norm(alfa-betaskh,"2")/norm(alfa,"2")
  MSEs[i] = MSEhk(X, h=1, k, sigma^2, beta, alfa)[[1]]
  NCs[i] = NCk(X,k)
  VIFs[i,] = VIFk(X,k)
  detRs[i] = detRk(X,k)
  i = i + 1
}
norms
MSEs
NCs
VIFs
detRs

##

minMSE(MSEs, discretizacion)

firstV(norms, threshold=0.1, discretizacion)
firstV(NCs, threshold=20, discretizacion)

firstV(VIFs[,1], threshold=10, discretizacion)
firstV(VIFs[,2], threshold=10, discretizacion)
firstVIF(VIFs, threshold=10, leap, stop)

firstV(detRs, threshold=0.1, discretizacion, minor = F)

###########################################################################
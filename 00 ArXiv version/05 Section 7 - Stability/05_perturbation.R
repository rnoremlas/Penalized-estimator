rm(list=ls())

source("03_functions.txt")
source("04_functions.txt")
source("05_functions.txt")

set.seed(2023)

############################################

datos = read.table("01_datos_multicol_alta.txt", header=T, sep=";")

Y = datos[,1]
X = as.matrix(datos[,-1])

##

reg = lm(Y~X+0) 
beta = as.double(reg$coefficients)
beta

####################################################################
# for a specific k 

h = 1
k = 1000

##
 
betahk = BetaKH(Y, X, h, k)[[2]]
betahk

betak = BetaKH(Y, X, h=0, k)[[2]]
betak

##

p = ncol(X)
alfa = matrix(mean(Y), p, 1)
for (i in 2:p) alfa[i] = cov(Y, X[,i])/var(X[,i])
BetaKH_bootstrap(Y, X, h, k, alfa)[[2]]

####################################################################
# stability

mu = 5
dv = 5
tol = 0.01

ite = 1000

norm_beta = numeric() 
norm_betak = numeric()
norm_betahk = numeric() 
n = nrow(X)
X.p = matrix(1, n, p)
for (i in 1:ite){
  for (j in 1:(p-1)) X.p[,j+1] = perturb(X[,j+1], mu, dv, tol) # keep in mind that in X there is a constant but not in X.p
  reg.p = lm(Y~X.p+0) 
  beta.p = as.double(reg.p$coefficients)
  betak.p = BetaKH_bootstrap(Y, X.p, h=0, k, alfa)[[2]]  # ridge estimator
  betahk.p = BetaKH_bootstrap(Y, X.p, h, k, alfa)[[2]]  # penalized estimator
  norm_beta[i] = (norm(beta-beta.p, "2")/norm(beta,"2"))*100
  norm_betak[i] = (norm(betak-betak.p, "2")/norm(betak,"2"))*100
  norm_betahk[i] = (norm(betahk-betahk.p, "2")/norm(betahk,"2"))*100
}

c(mean(norm_beta), quantile(norm_beta, 0.025), quantile(norm_beta, 0.975))
c(mean(norm_betak), quantile(norm_betak, 0.025), quantile(norm_betak, 0.975))
c(mean(norm_betahk), quantile(norm_betahk, 0.025), quantile(norm_betahk, 0.975))

####################################################################
# trace stability

leap = 1
discretization = seq(0, 200, leap)
norms = perturb_ite(Y, X, h, discretization, alfa, ite, mu, dv, tol)
write.table(norms, "05_norms_trace.txt", sep=";", col.names = T, row.names = F)

norm__beta = c(mean(norms[[2]]), quantile(norms[[2]], 0.025), quantile(norms[[2]], 0.975))
norm__ridge = c(mean(norms[[3]]), quantile(norms[[3]], 0.025), quantile(norms[[3]], 0.975))
norm__penalized = c(mean(norms[[4]]), quantile(norms[[4]], 0.025), quantile(norms[[4]], 0.975))
data.frame(norm__beta, norm__ridge, norm__penalized)

win.graph()
  plot(discretization, norms[[2]], type="l", col = "blue", lwd = 2, xlab="k", ylab="Stability of the OLS estimator")
  savePlot("05_Graph_StabilityOLS", type="eps")
  savePlot("05_Graph_StabilityOLS", type="pdf")
dev.off()

win.graph()
  plot(discretization, norms[[3]], type="l", col = "blue", lwd = 2, xlab="k", ylab="Stability of the ridge estimator")
  savePlot("05_Graph_StabilityRIDGE", type="eps")
  savePlot("05_Graph_StabilityRIDGE", type="pdf")
dev.off()

win.graph()
  plot(discretization, norms[[4]], type="l", col = "blue", lwd = 2, xlab="k", ylab="Stability of the penalized estimator")
  savePlot("05_Graph_StabilityPENALIZED", type="eps")
  savePlot("05_Graph_StabilityPENALIZED", type="pdf")
dev.off()

win.graph()
  plot(ts(cbind(norms[[2]], norms[[3]], norms[[4]]), start = 0, deltat = leap), plot.type="single", type="l", lty=c(1, 2, 3), col = c("orange", "lightblue", "blue"), lwd = 2, xlab="k", ylab="Stability of the OLS/ridge/penalized estimator")
  savePlot("05_Graph_StabilityRidgePENALIZED", type="eps")
  savePlot("05_Graph_StabilityRidgePENALIZED", type="pdf")
dev.off()


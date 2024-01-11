# en este código pretendo simular un ejemplo en el que si hay ortogonalidad, 
#   las estimaciones de la regresión multiple coincide con la de las simples.

rm(list = ls())

set.seed(2023)

#############################################################################
# multicolinealidad baja

obs = 100
cte = rep(1, obs)
X2 = rnorm(obs, 1, 10)
X3 = rnorm(obs, 1, 10)
u = rnorm(obs, 0, 2)

    X = cbind(cte, X2, X3)
    library(multiColl)
    CVs(X)
    VIF(X)
    CNs(X)
    RdetR(X)

beta1 = 5
beta2 = 2
beta3 = -4
betas = c(beta1, beta2, beta3)
Y = beta1 + beta2*X2 + beta3*X3 + u

datos_baja = cbind(Y, X)
write.table(datos_baja, "01_datos_multicol_baja.txt", row.names = F, sep=";")

reg.multiple = lm(Y~X2+X3)
#summary(reg.multiple)
betas.multiple = reg.multiple$coefficients

reg.simple2 = lm(Y~X2)
#summary(reg.simple2)
betas.simple2 = c(NA, reg.simple2$coefficients[[2]], NA)

reg.simple3 = lm(Y~X3)
#summary(reg.simple3)
betas.simple3 = c(NA, NA, reg.simple3$coefficients[[2]])

comparativa = data.frame(betas, betas.multiple, betas.simple2, betas.simple3)
comparativa


#############################################################################
# multicolinealida alta

for (i in 1: obs) X3[i] = 1 + 5*X2[i] + (-1)^i
  
X = cbind(cte, X2, X3)
library(multiColl)
CVs(X)
VIF(X)
CNs(X)
RdetR(X)

beta1 = 5
beta2 = 2
beta3 = -4
betas = c(beta1, beta2, beta3)
Y = beta1 + beta2*X2 + beta3*X3 + u

datos_alta = cbind(Y, X)
write.table(datos_alta, "01_datos_multicol_alta.txt", row.names = F, sep=";")

reg.multiple = lm(Y~X2+X3)
#summary(reg.multiple)
betas.multiple = reg.multiple$coefficients

reg.simple2 = lm(Y~X2)
#summary(reg.simple2)
betas.simple2 = c(NA, reg.simple2$coefficients[[2]], NA)

reg.simple3 = lm(Y~X3)
#summary(reg.simple3)
betas.simple3 = c(NA, NA, reg.simple3$coefficients[[2]])

comparativa = data.frame(betas, betas.multiple, betas.simple2, betas.simple3)
comparativa
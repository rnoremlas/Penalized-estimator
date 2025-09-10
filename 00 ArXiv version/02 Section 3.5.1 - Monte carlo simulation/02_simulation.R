rm(list=ls())

library(multiColl)
source("02_functions.txt")
library(dplyr)

############################################

set.seed(2023)

obs = seq(30,200,10)
  m1 = length(obs)
xis = c(0.96,0.97,0.98, 0.99)
  m2 = length(xis)
ps = 3:6
  m3 = length(ps)
des.tips = c(0.01, 0.1, 5, 10, 15) 
  m4 = length(des.tips)

m = m1*m2*m3*m4
m

############################################

res = matrix(, m, 13)
colnames(res) = c("n", "Xi", "p", "minCV", "maxVIF", "ECMkhmin", "ECMkh", "global_kh", "ECMkmin" , "ECMk", "global_k", "ECM", "CN")

top = 1
ite=0
for(n in obs){
  for(xi in xis){
    for(p in ps){
        for(des.tip in des.tips){
            ite = ite + 1
            #
            W = matrix(1, n, p)
            for (i in 2:p) W[,i] = rnorm(n, sample(c(-6, -4, -2, 0, 2 , 4, 6), 1), des.tip)
            #
            X = matrix(1 , n, p)
            for (i in 2:p) X[,i] = sqrt(1-xi^2)*W[,i] + xi*W[,p]
            #
            b = sample(c(-5, -4, -3, -2, -1, 1, 2, 3, 4, 5), p)
            u = rnorm(n, 0, 1)
            y = X%*%b+u 
            #
            reg = lm(y~X+0)
            beta = matrix(reg$coefficients, p, 1)
            sigma = summary(reg)[[6]]
            alfa = matrix(mean(y), p, 1)
            for (q in 2:p) alfa[q] = cov(y, X[,q])/var(X[,q])
            ECMkh = minMSEhk(X, h=1, sigma, beta, alfa, start = 0, leap = 0.01, stop = top, fig = 0)
            ECMk = minMSEhk(X, h=0, sigma, beta, alfa, start = 0, leap = 0.01, stop = top, fig = 0)
            ECM = MSEhk(X,h=0,k=0,sigma,beta,alfa)[[1]]
            #
            res[ite,1] = n
            res[ite,2] = xi
            res[ite,3] = p
            res[ite,4] = min(CVs(X))
            # if (p==2) res[ite,5] = 1
            # else res[ite,5] = max(VIF(X))
            res[ite,5] = max(VIF(X))
            res[ite,6] = ECMkh[[1]]
            res[ite,7] = ECMkh[[2]]
            res[ite,8] = ECMkh[[3]]
            res[ite,9] = ECMk[[1]]
            res[ite,10] = ECMk[[2]]
            res[ite,11] = ECMk[[3]]
            res[ite,12] = ECM
            res[ite,13] = CN(X)
            #
            print(ite)
        }
    }
  }
}
ite
#res

write.table(res, file=paste("02_datos_simulados_top=", top,".txt", sep=""), row.names=F, sep=";") 
# res = read.table("02_datos_simulados.txt", header=F)

#######################################################

minCV = res[,4]
  summary(minCV)
  
maxVIF = res[,5]
  summary(maxVIF)
  
CN = res[,13]
  summary(CN)  

n = res[,1]
p = res[,3]    
ECMkh = res[,7]
ECMk = res[,9]
ECM = res[,12]
global_kh = res[,8]
global_h = res[,11]
output = data.frame(ECMkh, ECMk, ECM, minCV, maxVIF, CN, global_kh, global_h, n, p)
total1 = nrow(output)

# 

output1 = filter(output, (global_kh >0) & (global_h > 0))
  write.table(output1, file=paste("02_datos_simuladoss_filtrados_top=", top,".txt", sep=""), row.names=F, sep=";")
  summary(output1)
  total2 = nrow(output1)
  (total2/total1)*100 # percentage of simulated cases where the minimum ECM can be calculated in penalized and ridge estimator
classification(output1)

outputA = filter(output1, (ECMkh < ECMk) & (ECMk < ECM))
  nrow(outputA)
  summary(outputA)

outputB = filter(output1, (ECMkh < ECM) & (ECM < ECMk))
  nrow(outputB)
  summary(outputB)

outputC = filter(output1, (ECMk < ECMkh) & (ECMkh < ECM))
  nrow(outputC)
  summary(outputC)
  
##########################################################

install.packages("glmnet")
install.packages("ggplot2")
install.packages("ElemStatLearn")
install.packages("parcor")
install.packages("pracma")
install.packages("PreProcess")
install.packages("Rtools")
install.packages("oompaBase")
install.packages("glmnet")
library(oompaBase)
library(PreProcess)
library(multiColl)
library(ElemStatLearn)
library(parcor)
library(multiColl)
library(glmnet)
require(glmnet)
require(ggplot2)
require(ElemStatLearn)
require(parcor)
require(pracma)

#funciones del penalizado
source("02_functions.txt")
source("03_functions.txt")
source("04_functions.txt")
source("05_functions.txt")
source("06_functions.txt")


# Cargar datos libreriaElemStatLearn
data=data(prostate)
attach(prostate)


# Construcción de matrices
n <- nrow(prostate)
cte = rep(1,n)
X = as.matrix(cbind(prostate[,-((ncol(prostate)-1):ncol(prostate))]))
X_o= cbind(cte, X)
y <- prostate$lpsa
p = ncol(X_o)
CVs(X_o)
VIF(X_o)
CNs(X_o)
RdetR(X_o)

# Inicializar vectores para almacenar los RMSE
num_iter <- 100
rss_results <- matrix(NA, nrow = num_iter, ncol = 4 * 2)#4 metodos train y test

set.seed(123) 
for (i in 1:100) {
  # Partición 70%-30%
  train_idx <- sample(1:n, size = floor(0.7 * n))
  test_idx <- setdiff(1:n, train_idx)
  
  X_train <- X_o[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X_o[test_idx, ]
  y_test <- y[test_idx]


##APLICACION OLS (DATOS NO ESTANDARIZADOS)
reg <- lm(y_train ~ X_train - 1)
beta <- matrix(reg$coef, p, 1)
sigma <- as.numeric(summary(reg)[6])

## OBTENCION DE ALFA
alfa <- matrix(mean(y_train), p, 1)
for (j in 2:p) {
  alfa[j] <- cov(y_train, X_train[, j]) / var(X_train[, j])
}
# para comprobar hacer esto: lm(y_train~X_train[,3])

##############penalizado#####################
# OBTENCION DE TRAZA (PENALIZADO)

start = 0
leap = 0.01
stop = 5

h=1 # h=0 cresta con h=1 penalizado

##

output = TracesHK(y_train, X_train, alfa, h, sigma, beta, start, leap, stop, graph=0)
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
perturb_ite(y_train, X_train, h, discretization=0.01, alfa, ite, mu=5, dv=5, tol=0.01)

#####################################################
# choice of k (SOLO PARA EL MIN MSE)

firstV(norms.rate, threshold=0.1, kas)

KminMSE=as.numeric(minMSE(MSEs, kas)[1])

#####################################################

# ESTIMACION for a fixed k

P=PenalizedEstimation(y_train, X_train, h, KminMSE, alfa, sigma, beta, ite, tol = 0.01, setseed = 1)
beta_pen=as.numeric(P[1:p,2])

y_hat_train_pen = X_train %*% beta_pen
rss_results[i, 1] = sqrt(sum((y_train - y_hat_train_pen)^2)/n)
y_hat_test_pen <- X_test %*% beta_pen
rss_results[i, 2] = sqrt(sum((y_test - y_hat_test_pen)^2)/n)

################Cresta (con paquete penalizado)##############
# OBTENCION DE TRAZA (ridge)

start = 0
leap = 0.01
stop = 5

h=0 # h=0 cresta con h=1 penalizado

##

output = TracesHK(y_train, X_train, alfa, h, sigma, beta, start, leap, stop, graph=0)
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
perturb_ite(y_train, X_train, h, discretization=0.01, alfa, ite, mu=5, dv=5, tol=0.01)

#####################################################
# choice of k (SOLO PARA EL MIN MSE)

firstV(norms.rate, threshold=0.1, kas)

KminMSE=as.numeric(minMSE(MSEs, kas)[1])

#####################################################

# ESTIMACION for a fixed k

P=PenalizedEstimation(y_train, X_train, h, KminMSE, alfa, sigma, beta, ite, tol = 0.01, setseed = 1)
beta_ridge=as.numeric(P[1:p,2])
y_hat_train_ridge = X_train %*% beta_ridge
rss_results[i, 3]= sqrt(sum((y_train - y_hat_train_ridge)^2)/n)

y_hat_test_ridge <- X_test %*% beta_ridge
rss_results[i, 4] = sqrt(sum((y_test - y_hat_test_ridge)^2)/n)


########################ESTIMACION lasso MEDIANTE GLMNET###########

lasso_model <- cv.glmnet(X_train, y_train, alpha = 1, grouped = FALSE, standardize = FALSE)
best_lambda_lasso <- lasso_model$lambda.min
y_lasso_pred_train <- predict(lasso_model, X_train, s = best_lambda_lasso)
rss_results[i, 5]=sqrt(sum((y_train - y_lasso_pred_train)^2) / n)

y_lasso_pred_test <- predict(lasso_model, X_test, s = best_lambda_lasso)
rss_results[i, 6]=sqrt(sum((y_test - y_lasso_pred_test)^2) / n)


########################ESTIMACION elastic net MEDIANTE GLMNET###########

elasticnet_model <- cv.glmnet(X_train, y_train, alpha = 0.5, grouped = FALSE, standardize = FALSE)
best_lambda_elasticnet <- elasticnet_model$lambda.min
y_elasticnet_pred_train <- predict(elasticnet_model, X_train, s = best_lambda_elasticnet)
rss_results[i, 7]=sqrt(sum((y_train - y_elasticnet_pred_train)^2) / n)

y_elasticnet_pred_test <- predict(elasticnet_model, X_test, s = best_lambda_elasticnet)
rss_results[i, 8]=sqrt(sum((y_test - y_elasticnet_pred_test)^2) / n)
}

# MEAN
rss_means <- colMeans(rss_results, na.rm = TRUE)
rss_table_prostate <- data.frame(
  Method = c("GRR", "Ridge", "Lasso", "Elastic Net"),
  RSS_Train = rss_means[c(1, 3, 5, 7)],
  RSS_Test = rss_means[c(2, 4, 6, 8)]
)

print(rss_table_prostate)
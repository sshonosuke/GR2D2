###-------------------------------------------------------------###
###             Model validation for abalone data               ###
###-------------------------------------------------------------###
rm(list=ls())

## R packages
library(fda)                  # for basis functions
source("gR2D2.R")             # grouped R2D2 prior 
source("Competitors.R")       # other shrinkage priors

## Dataset 
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/abalone/abalone.data"
abalone <- read.csv(url, header=F)
colnames(abalone) <- c("Sex", "Length", "Diameter", "Height", "WholeWeight", 
                       "ShuckedWeight", "VisceraWeight", "ShellWeight", "Rings")

data_M <- subset(abalone, Sex=="M")
Y <- data_M$Rings
X_origin <- cbind(data_M[,2:8])
G <- dim(X_origin)[2]      # number of groups (= number of additive effects)
n <- dim(X_origin)[1]      # sample size

## basis function (B spline) 
L <- 10        # number of basis functions
ID <- rep(1:G, rep(L, G))
p <- length(ID)
X <- matrix(NA, n, p)

Basis <- list()
for(g in 1:G){
  basis_obj <- create.bspline.basis(rangeval=range(X_origin[,g]), nbasis=L)
  X[,ID==g] <- eval.basis(X_origin[,g], basis_obj)
  Basis[[g]] <- basis_obj
}


## out-of-sample validation
set.seed(1)

m <- 500   # number of test samples
R <- 100   # number of Monte Carlo replications

Method <- c("R2D2", "HS", "GIGG", "GHS", "gR2D2-D", "gR2D2-L")
M <- length(Method)
MSE <- CP <- AL <- matrix(NA, R, M)

# Replications 
for(r in 1:R){
  try({
    test_id <- sort(sample(1:n, m))
    X_test <- X[test_id,]
    X_train <- X[-test_id,]
    Y_test <- Y[test_id]
    Y_train <- Y[-test_id]
    
    # MCMC
    Fit <- list()
    mc <- 3000
    bn <- 1000
    Fit[[1]] <- R2D2(Y_train, X_train, a=1/2, b=1/2, mc=mc, bn=bn) 
    Fit[[2]] <- HS(Y_train, X_train, a=1/2, b=1/2, mc=mc, bn=bn) 
    Fit[[3]] <- GIGG(Y_train, X_train, ID, a=1/p, b=1/2, mc=mc, bn=bn) 
    Fit[[4]] <- GHS(Y_train, X_train, ID, mc=mc, bn=bn) 
    Fit[[5]] <- R2D2_group(Y_train, X_train, Group_ID=ID, phi_dist="Dir", mc=mc, bn=bn) 
    Fit[[6]] <- R2D2_group(Y_train, X_train, Group_ID=ID, phi_dist="LN", mc=mc, bn=bn) 
    
    # prediction 
    Pred <- matrix(NA, m, M)
    PI <- array(NA, c(2, m, M))
    for(k in 1:M){
      ep <- Fit[[k]]$sig * matrix(rnorm(mc*m), mc, m)
      rn <- X_test%*%t(Fit[[k]]$beta) + t(ep)
      Pred[,k] <- apply(rn, 1, mean)
      PI[,,k] <- apply(rn, 1, quantile, prob=c(0.025, 0.975))
    }
    
    MSE[r,] <- apply((Pred-Y_test)^2, 2, mean)
    AL[r,] <- apply(PI[2,,]-PI[1,,], 2, mean)
    CP[r,] <- apply(PI[2,,]>Y_test & PI[1,,]<Y_test, 2, mean)
    
    print(r)
  })
}


## Table 3
mse <- apply(na.omit(MSE), 2, mean)[3:6]
cp <- 100*apply(na.omit(CP), 2, mean)[3:6]
tab <- cbind(t(AL_reg), mse, cp)
write.csv(tab, file="Table3.csv")


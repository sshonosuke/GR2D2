###-------------------------------------------------------------###
###             Example of grouped R2D2 prior                   ###
###-------------------------------------------------------------###
## Preparation
rm(list=ls())
source("gR2D2.R")
source("Competitors.R")
library(gigg)
library(MASS)
mat_fun <- function(p, ID, rin=0.7, rout=0.2){
  mat <- matrix(rout, p, p)
  for(k in 1:max(ID)){
    mat[ID==ID[k], ID==ID[k]] <- rin
  }
  diag(mat) <- 1
  return(mat)
}

meth <- c("R2D2", "HS", "GIGG", "GHS", "GR2D2-D", "GR2D2-L")



## settings
scenario <- 1   # 1-5
n <- 250
p <- 50
SNR <- 0.7

# true beta
if(scenario==1){
  Beta <- rep(0, p)
  G <- 5
  p_G <- rep(10, G)
  ID <- rep(1:G, p_G)
  Beta[c(1,11,21,31,41)] <- c(5, -5, 5, -5, 5)
  mat <- mat_fun(p, ID)
  signal <- (1:p)[Beta!=0]
  sum(abs(Beta))
}

if(scenario==2){
  Beta <- rep(0, p)
  G <- 5
  p_G <- rep(10, G)
  ID <- rep(1:G, p_G)
  Beta[1:3] <- 3
  Beta[11:13] <- (-3)
  Beta[21] <- 2
  Beta[41:45] <- (-1)
  mat <- mat_fun(p, ID)
  signal <- (1:p)[Beta!=0]
  sum(abs(Beta))
}


if(scenario==3){
  Beta <- rep(0, p)
  G <- 5
  p_G <- rep(10, G)
  ID <- rep(1:G, p_G)
  Beta[1:3] <- 3
  Beta[11:13] <- (-3)
  Beta[c(41,45)] <- c(3,4)
  mat <- mat_fun(p, ID)
  signal <- (1:p)[Beta!=0]
  sum(abs(Beta))
}


if(scenario==4){
  Beta <- rep(0, p)
  G <- 5
  p_G <- rep(10, G)
  ID <- rep(1:G, p_G)
  Beta[c(1,3,5,7,9)] <- 3
  Beta[21:25] <- (-2)
  mat <- mat_fun(p, ID)
  signal <- (1:p)[Beta!=0]
  sum(abs(Beta))
}

if(scenario==5){
  Beta <- rep(0, p)
  G <- 5
  p_G <- rep(10, G)
  ID <- rep(1:G, p_G)
  Beta[c(1,3,5,7,9)] <- c(5, -5, 5, -5, 5)
  mat <- mat_fun(p, ID)
  signal <- (1:p)[Beta!=0]
  sum(abs(Beta))
}


## data generation
# covariate 
X <- mvrnorm(n, rep(0,p), mat)
Reg_variation <- as.vector( t(Beta)%*%mat%*%Beta )

# error variance 
Sig <- sqrt( (1-SNR)/SNR*Reg_variation )

# response
Y <- as.vector(X%*%Beta) + Sig*rnorm(n)
Y <- Y - mean(Y)


## MCMC
mc <- 3000
bn <- 1000
Fit <- list()
Fit[[1]] <- R2D2(Y, X, a=1/2, b=1/2, mc=mc, bn=bn)       # standard R2D2 prior
Fit[[2]] <- HS(Y, X, a=1/2, b=1/2, mc=mc, bn=bn)         # standard horseshoe prior
Fit[[3]] <- GIGG(Y, X, ID, a=1/p, b=1/2, mc=mc, bn=bn)   # GIGG prior
Fit[[4]] <- GHS(Y, X, ID, mc=mc, bn=bn)                  # grouped horseshoe prior 


variation <- c()
for(g in 1:G){
  variation[g] <- mean( predict(lm(Y~X[,ID==g]))^2 )
}
ratio <- variation/sum(variation)
c <- 1/2
Fit[[5]] <- R2D4_Laplace(Y, X, Group_ID=ID, a_g=c*ratio, b=c, phi_dist="Dir", mc=mc, bn=bn)   # grouped R2D2 (Dirichlet)
Fit[[6]] <- R2D4_Laplace(Y, X, Group_ID=ID, a_g=c*ratio, b=c, phi_dist="LN", mc=mc, bn=bn)    # grouped R2D2 (Logistic normal)

Fit[[6]]$ac_ratio      # acceptance ratio 


## point estimate
L <- length(Fit)
hBeta <- c()
for(k in 1:L){
  hBeta <- cbind(hBeta, apply(Fit[[k]]$beta, 2, median))
}
dimnames(hBeta)[[2]] <- meth


## credible interval 
CI <- list()
for(k in 1:L){
  CI[[k]] <- apply(Fit[[k]]$beta, 2, quantile, prob=c(0.025, 0.975))
}

## coverage and average length (overall)
CP_overall <- AL <- c()
for(k in 1:L){
  CP_overall[k] <- mean(CI[[k]][1,]<Beta & Beta<CI[[k]][2,]) 
  AL[k] <- mean(CI[[k]][2,]-CI[[k]][1,])
}
names(CP_overall) <- names(AL) <- meth
100*CP_overall
AL

# MSE (total)
apply((hBeta-Beta)^2, 2, sum)


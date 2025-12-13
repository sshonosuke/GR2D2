###-------------------------------------------------------------###
###                Monte Carlo simulation                       ###
###      of Yanchenko, Irie and Sugasawa (2024, arXiv)          ###
###-------------------------------------------------------------###
rm(list=ls())

## settings
set.seed(1)
scenario <- 3    # scenario of true regression coefficients (from 1 to 5)
p <- 200          # number of covariates (50 / 200 / 1000)

# sample size (dependent on 'p')
if(p==50){ n <- 250 }
if(p==200){ n <- 150 }
if(p==1000){ n <- 200 }

# number of Monte Carlo replications
if(p==50 | p==200){ R <- 500 }    
if(p==1000){ R <- 100 }    

# group structure
G <- p/10
p_G <- rep(p/G, G)
ID <- rep(1:G, p_G)


## load packages
source("gR2D2.R")
source("Competitors.R")
library(MASS)

## function for correlation matrix
mat_fun <- function(p, ID, rin=0.7, rout=0.2){
  mat <- matrix(rout, p, p)
  for(k in 1:max(ID)){
    mat[ID==ID[k], ID==ID[k]] <- rin
  }
  diag(mat) <- 1
  return(mat)
}

## True regression coefficients
if(scenario==1){
  Beta <- rep(0, p)
  Beta[c(1,11,21,31,41)] <- c(5, -5, 5, -5, 5)
  mat <- mat_fun(p, ID)
  signal <- (1:p)[Beta!=0]
}

if(scenario==2){
  Beta <- rep(0, p)
  Beta[1:3] <- 3
  Beta[11:13] <- (-3)
  Beta[21] <- 2
  Beta[41:45] <- (-1)
  mat <- mat_fun(p, ID)
  signal <- (1:p)[Beta!=0]
}

if(scenario==3){
  Beta <- rep(0, p)
  Beta[1:3] <- 3
  Beta[11:13] <- (-3)
  Beta[c(41,45)] <- c(3,4)
  mat <- mat_fun(p, ID)
  signal <- (1:p)[Beta!=0]
}

if(scenario==4){
  Beta <- rep(0, p)
  Beta[c(1,3,5,7,9)] <- 3
  Beta[21:25] <- (-2)
  mat <- mat_fun(p, ID)
  signal <- (1:p)[Beta!=0]
}

if(scenario==5){
  Beta <- rep(0, p)
  Beta[c(1,3,5,7,9)] <- c(5, -5, 5, -5, 5)
  mat <- mat_fun(p, ID)
  signal <- (1:p)[Beta!=0]
}


## Preparation for Monte Carlo replications
SNR_set <- c(5,7)/10      # two scenarios for SNR
L <- length(SNR_set)
Method <- c("R2D2", "HS", "GIGG", "GHS", "gR2D2-D", "gR2D2-LN")
M <- length(Method)

MSE_signal <- MSE_noise <- array(NA, c(R, M, L))
CP_all <- AL_all <- array(NA, c(R, M, L))
AC <- matrix(NA, R, L)


## Main part (Monte Carlo replication)
mc <- 3000     # number of posterior samples after burn-in 
bn <- 1000     # burn-in length 

for(l in 1:L){
  SNR <- SNR_set[l]
  
  # sigma (depends on SNR)
  Reg_variation <- as.vector( t(Beta)%*%mat%*%Beta )
  Sig <- sqrt( (1-SNR)/SNR*Reg_variation )
  
  # replications
  Fit <- list()
  for(r in 1:R){
    # covariate
    X <- mvrnorm(n, rep(0, p), mat)
    
    # response variables 
    Y <- as.vector(X%*%Beta) + Sig*rnorm(n)
    Y <- Y - mean(Y)
    
    # MCMC
    try( Fit[[1]] <- R2D2(Y, X, a=1/2, b=1/2, mc=mc, bn=bn) )       # R2D2 (without group information)
    try( Fit[[2]] <- HS(Y, X, a=1/2, b=1/2, mc=mc, bn=bn) )         # HS (without group information)
    try( Fit[[3]] <- GIGG(Y, X, ID, a=1/p, b=1/2, mc=mc, bn=bn) )   # grouped GIG
    try( Fit[[4]] <- GHS(Y, X, ID, mc=mc, bn=bn) )                  # grouped HS 
    try( Fit[[5]] <- R2D2_group(Y, X, Group_ID=ID, phi_dist="Dir", mc=mc, bn=bn) )   # (proposed) grouped R2D2 with Dirichlet 
    try( Fit[[6]] <- R2D2_group(Y, X, Group_ID=ID, phi_dist="LN", mc=mc, bn=bn) )    # (proposed) grouped R2D2 with logistic normal
    
    # posterior median
    hBeta <- c()
    for(k in 1:M){
      hBeta <- cbind(hBeta, apply(Fit[[k]]$beta, 2, median))
    }
    
    # MSE (signal)
    MSE_signal[r,,l] <- apply((hBeta[signal,]-Beta[signal])^2, 2, sum)
    
    # MSE (noise)
    MSE_noise[r,,l] <- apply((hBeta[-signal,]-Beta[-signal])^2, 2, sum)
    
    # acceptance rate 
    AC[r,l] <- mean( Fit[[6]]$ac_ratio )
    
    # 95% credible interval 
    CI <- list()
    for(k in 1:M){
      CI[[k]] <- apply(Fit[[k]]$beta, 2, quantile, prob=c(0.025, 0.975))
    }
    
    # coverage and length
    for(k in 1:M){
      CP_all[r,k,l] <- mean(CI[[k]][1,]<Beta & Beta<CI[[k]][2,]) 
      AL_all[r,k,l] <- mean(CI[[k]][2,]-CI[[k]][1,]) 
    }
    print(r)
  }
  
  # print
  print(l)
}


## save
save(list=ls(), file=paste0("Sim-p", p, "-result", scenario, ".RData"))



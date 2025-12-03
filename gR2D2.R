library(MCMCpack)
library(GIGrvg)
library(statmod)
library(base)   # for trigamma function 


###---------------------------------------------------------------###
###           Grouoed R2D2 prior with Laplace kernel              ###  
###---------------------------------------------------------------###   
## phi_dist: distribution of phi 
## ("Dir" for Dirichlet distribution and "LN" for logistic normal) 
## eta_sig: standard deviation of normal distributions for LN

R2D4_Laplace <- function(Y, X, Group_ID, mc=2000, bn=1000, a_g=NULL, b=1/2, 
                         phi_dist="LN"){
  # preparation
  p <- dim(X)[2]
  G <- max(Group_ID)
  X_mat <- t(X)%*%X
  n0 <- d0 <- 2        # prior of sigma^2: IG(n0/2, d0/2)
  thres <- 10^(-10)
  p_g <- table(Group_ID)    # number of covariates in each group
  
  # function: log-ratio-transformation 
  LR <- function(x){
    K <- length(x)
    log(x[1:(K-1)]/x[K])
  }
  
  # tuning parameters
  sID <- rep(1:G, p_g-1)    # ID for eta (only for "LN")
  if(is.null(a_g)){ a_g <- rep(1/2, G) } 
  a_gj <- a_g[Group_ID]/p_g[Group_ID]     # tuning parameter for "Dir" and "LN"
  if(phi_dist=="LN"){
    IS_eta <- list()
    for(g in 1:G){
      mat <- solve( trigamma(1/p_g[g])*matrix(1, p_g[g], p_g[g]) + trigamma(1/p_g[g])*diag(p_g[g]) )
      IS_eta[[g]] <- mat[-p_g[g], -p_g[g]]
    }
  }
  
  # initial values
  beta <- coef(lm(Y~X-1))
  sig  <- sd(Y-X%*%beta)
  psi <- rep(1, p)  # within-group parameter for Laplace distribution 
  phi <- rep(1/p_g, p_g)  # within-group parameter for local shrinkage
  gam <- rep(1, G)  # group-wise parameter for group-level shrinkage
  tau <- 1          # global variance for shrinkage 
  
  # matrix to store posterior samples
  beta_pos <- matrix(NA, mc, p)
  sig_pos <- rep(NA, mc)
  tau_pos <- rep(NA, mc)
  gam_pos <- matrix(NA, mc, G)
  phi_pos <- matrix(NA, mc, p)
  AC <- matrix(0, mc, G)    # acceptance indicator (for LN)
  
  # MCMC
  MC <- mc + bn
  for(k in 1:MC){
    # beta
    lam <- psi*phi*gam[Group_ID]/2
    lam[lam<thres] <- thres
    Inv_Lam <- diag( 1/lam )
    Sig_beta <- solve(  X_mat + Inv_Lam )
    Mu_beta <- Sig_beta%*%t(X)%*%Y
    beta <- mvrnorm(1, Mu_beta, sig^2*Sig_beta)
    
    # sigma (error standard deviation)
    resid <- Y - as.vector(X%*%beta)
    n1 <- n0 + n + p
    d1 <- d0 + sum(resid^2) + sum(beta^2/lam)
    sig <- sqrt( rinvgamma(1, n1/2, d1/2) )
    
    # psi (within-group parameter for Laplace)
    psi_mu <- sqrt( sig^2*phi*gam[Group_ID]/2 ) / abs(beta)
    psi <- 1/rinvgauss(p, psi_mu, 1)
    
    ## local parameters (Dir)
    if(phi_dist=="Dir"){
      # phi (within-group parameter for local shrinkage)
      T_gj <- c()
      for(j in 1:p){
        chi_g <- beta[j]^2 / ( sig^2*psi[j]*gam[Group_ID[j]]/2 )
        T_gj[j] <- rgig(1, lambda=a_gj[j]-0.5, chi=chi_g, psi=2)
      }
      for(g in 1:G){
        phi[Group_ID==g] <- T_gj[Group_ID==g] / sum(T_gj[Group_ID==g]) 
      }
      # gamma (group-level local parameters)
      for(g in 1:G){
        beta_g <- beta[Group_ID==g]
        psi_g <- psi[Group_ID==g]
        phi_g <- phi[Group_ID==g]
        chi_gam <- sum( beta_g^2 / (sig^2*psi_g*phi_g/2) ) 
        lam_gam <- a_g[g] - p_g[g]/2 
        gam[g] <- rgig(1, lambda=lam_gam, chi=chi_gam, psi=2/tau^2) 
      }
    }
    
    ## local parameters (LN)
    if(phi_dist=="LN"){
      # phi & gamma (group-level and within-group parameter)
      T_gj <- c()
      for(j in 1:p){
        chi_g <- beta[j]^2 / ( sig^2*psi[j]*gam[Group_ID[j]]/2 )
        T_gj[j] <- rgig(1, lambda=a_gj[j]-0.5, chi=chi_g, psi=2)
      }
      for(g in 1:G){
        # phi
        phi_new <- T_gj[Group_ID==g] / sum(T_gj[Group_ID==g])   # proposal of phi
        # gamma
        beta_g <- beta[Group_ID==g]
        psi_g <- psi[Group_ID==g]
        phi_g <- phi[Group_ID==g]
        chi_gam <- sum( beta_g^2 / (sig^2*psi_g*phi_g/2) ) 
        lam_gam <- a_g[g] - p_g[g]/2 
        gam_new <- rgig(1, lambda=lam_gam, chi=chi_gam, psi=2/tau^2)   # proposal of gamma
        # accept/reject
        phi_sub <- phi[Group_ID==g]
        log_f_old <- (-1)*sum( a_gj[Group_ID==g]*log(phi_sub) ) - (0.5)*t(LR(phi_sub))%*%IS_eta[[g]]%*%LR(phi_sub)
        log_f_new <- (-1)*sum( a_gj[Group_ID==g]*log(phi_new) ) - (0.5)*t(LR(phi_new))%*%IS_eta[[g]]%*%LR(phi_new)
        if(log(runif(1))<(log_f_new-log_f_old)){
          phi[Group_ID==g] <- phi_new
          gam[g] <- gam_new
          if(k > bn){ AC[k-bn, g] <- 1 }
        }
      }
    }
    
    ## local parameters (fixed phi)
    if(phi_dist=="fix"){
      # gamma (group-level local parameters)
      for(g in 1:G){
        beta_g <- beta[Group_ID==g]
        psi_g <- psi[Group_ID==g]
        phi_g <- phi[Group_ID==g]
        chi_gam <- sum( beta_g^2 / (sig^2*psi_g*phi_g/2) ) 
        lam_gam <- a_g[g] - p_g[g]/2 
        gam[g] <- rgig(1, lambda=lam_gam, chi=chi_gam, psi=2/tau^2) 
      }
    }
    
    # tau (global shrinkage parameter)
    a_tau <- b + sum(a_g)
    b_tau <- 1 + sum(gam) 
    tau <- sqrt( rinvgamma(1, a_tau, b_tau) )
    
    # store MCMC samples
    if(k > bn){
      beta_pos[k-bn,] <- beta
      sig_pos[k-bn] <- sig
      tau_pos[k-bn] <- tau
      gam_pos[k-bn,] <- gam
      phi_pos[k-bn,] <- phi
    }
  }
  
  # output
  AC_ratio <- apply(AC, 2, mean)
  Result <- list(beta=beta_pos, sig=sig_pos, tau=tau_pos, phi=phi_pos, 
                 gam=gam_pos, ac_ratio=AC_ratio)
  return(Result)
}





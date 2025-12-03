###-------------------------------------------------------------###
###             Comparative shrinkage priors                    ###
###-------------------------------------------------------------###

## standard Horseshoe 
HS <- function(Y, X, a=1/2, b=1/2, mc=2000, bn=500){
  # preparation
  p <- dim(X)[2]
  X_mat <- t(X)%*%X
  delta <- 1   # prior of sigma^2
  thres <- 10^(-8)
  
  # initial values
  beta <- coef(lm(Y~X-1))
  sig <- 1
  psi <- rep(1, p)
  xi <- 1
  
  # matrix to store posterior samples
  beta.pos <- matrix(NA, mc, p)
  sig.pos <- c()
  
  # MCMC
  for(k in 1:(mc+bn)){
    # beta
    ss <- psi*xi
    ss[ss<thres] <- thres
    A <- solve( X_mat + diag(1/ss) )
    B <- as.vector(t(X)%*%Y)
    beta <- mvrnorm(1, A%*%B, sig^2*A)
    if(k > bn){
      beta.pos[k-bn,] <- beta
    }
    
    # sig
    resid <- Y - as.vector(X%*%beta)
    sig <- sqrt( rinvgamma(1, delta+(n+p)/2, delta+0.5*sum(resid^2)+0.5*sum(beta^2/(psi*xi)) ) )
    
    if(k > bn){
      sig.pos[k-bn] <- sig
    }
    
    # psi (local parameters)
    nu <- rinvgamma(p, 1, 1+1/psi)
    psi <- rinvgamma(p, 1, 1/nu + beta^2/(2*sig^2*xi))
    psi[psi<thres] <- thres
    
    # xi (global parameter)
    gam <- rinvgamma(1, 1, 1+1/xi)
    xi <- rinvgamma(1, (p+1)/2, 1/gam + sum(beta^2/(2*sig^2*psi)))
    xi[xi<thres] <- thres
  }
  return(list(beta=beta.pos, sig=sig.pos))
}






## standard R2D2
R2D2 <- function(Y, X, a=1/2, b=1/2, mc=2000, bn=500){
  # preparation
  p <- dim(X)[2]
  X_mat <- t(X)%*%X
  delta <- 1   # prior of sigma^2
  thres <- 10^(-8)
  
  # initial values
  beta <- coef(lm(Y~X-1))
  sig <- 1
  lam <- rep(1, p)
  psi <- rep(1, p)
  xi <- 1
  api <- a/p
  
  # matrix to store posterior samples
  beta.pos <- matrix(NA, mc, p)
  sig.pos <- c()
  
  # MCMC
  for(k in 1:(mc+bn)){
    # beta
    A <- solve( X_mat + diag(2/(psi*lam)) )
    B <- as.vector(t(X)%*%Y)
    beta <- mvrnorm(1, A%*%B, sig^2*A)
    if(k > bn){
      beta.pos[k-bn,] <- beta
    }
    
    # sig
    resid <- Y - as.vector(X%*%beta)
    sig <- sqrt( rinvgamma(1, delta+(n+p)/2, delta+0.5*sum(resid^2)+0.5*sum(beta^2/(psi*lam/2)) ) )
    
    if(k > bn){
      sig.pos[k-bn] <- sig
    }
    
    # psi
    for(j in 1:p){
      psi[j] <- 1/rinvgauss(1, sqrt(sig^2*lam[j]/2)/abs(beta[j]))
    }
    
    psi[psi<thres] <- thres
    psi <- 1/psi
    
    # local parameters
    for(j in 1:p){
      lam[j] <- rgig(1, lambda=api-1/2, chi=2*beta[j]^2/(psi[j]*sig^2), psi=2*xi)
    }
    lam[lam<thres] <- thres
    
    xi <- rgamma(1, a+b, 1+sum(lam))
  }
  return(list(beta=beta.pos, sig=sig.pos))
}




## GIGG
GIGG <- function(Y, X, ID, a=1/2, b=1/2, mc=2000, bn=500){
  # preparation
  G <- max(ID)
  p <- dim(X)[2]
  X_mat <- t(X)%*%X
  delta <- 1   # prior of sigma^2
  thres <- 10^(-10)
  
  # initial values
  beta <- coef(lm(Y~X-1))
  sig <- 1
  tau <- 1
  lam <- rep(1, p)
  gam <- rep(1, G)
  
  # matrix to store posterior samples
  beta.pos <- matrix(NA, mc, p)
  sig.pos <- c()
  
  # MCMC
  for(k in 1:(mc+bn)){
    # beta
    V <- tau^2*lam*gam[ID]
    V[V<thres] <- thres
    A <- solve( X_mat/sig^2 + diag(1/V) )
    B <- as.vector(t(X)%*%Y)/sig^2
    beta <- mvrnorm(1, A%*%B, A)
    if(k > bn){
      beta.pos[k-bn,] <- beta
    }
    
    # sig
    resid <- Y - as.vector(X%*%beta)
    sig <- sqrt( rinvgamma(1, delta+n/2, delta+0.5*sum(resid^2)) )
    if(k > bn){
      sig.pos[k-bn] <- sig
    }
    
    # lambda
    lam <- rinvgamma(p, b+1/2, 1+beta^2/(2*tau^2*gam[ID]))
    
    # gam
    for(g in 1:G){
      gam[g] <- rgig(1, lambda=a-1/2, chi=sum(beta[ID==g]^2/(tau^2*lam[ID==g])), psi=2)
    }
    
    # tau
    tau <- sqrt( rinvgamma(1, delta+0.5*p, delta+0.5*sum(beta^2/(gam[ID]*lam))) )
  }
  return(list(beta=beta.pos, sig=sig.pos))
}




## GHS
GHS <- function(Y, X, ID, mc=2000, bn=500){
  # preparation
  G <- max(ID)
  p <- dim(X)[2]
  X_mat <- t(X)%*%X
  delta <- 1   # prior of sigma^2
  thres <- 10^(-10)
  
  # initial values
  beta <- coef(lm(Y~X-1))
  sig <- 1
  tau <- 1
  lam <- rep(1, G)
  tt <- rep(1, G)
  delta <- rep(1, p)
  cc <- rep(1, p)
  
  # matrix to store posterior samples
  beta.pos <- matrix(NA, mc, p)
  sig.pos <- c()
  
  # MCMC
  for(k in 1:(mc+bn)){
    # beta
    V <- tau^2*lam[ID]*delta
    A <- solve( X_mat + diag(1/V) + 10^(-5)*diag(p) )
    B <- as.vector(t(X)%*%Y)
    beta <- mvrnorm(1, A%*%B, sig^2*A)
    if(k > bn){
      beta.pos[k-bn,] <- beta
    }
    
    # sig
    resid <- Y - as.vector(X%*%beta)
    sig <- sqrt( rinvgamma(1, (n-1+p)/2, 0.5*sum(resid^2)+0.5*sum(beta^2/V)) )
    if(k > bn){
      sig.pos[k-bn] <- sig
    }
    
    # lambda
    for(g in 1:G){
      sg <- sum(ID==g)
      tt[g] <- rinvgamma(1, 1, 1+1/lam[g])
      lam[g] <- rinvgamma(1, (sg+1)/2, 0.5*sum(beta[ID==g]^2/(sig^2*tau^2*delta[ID==g]))+1/tt[g])
    }
    
    # delta 
    cc <- rinvgamma(1, 1, 1+1/delta)
    delta <- rinvgamma(p, 1, 0.5*beta^2/(sig^2*tau^2*lam[ID])+1/cc)
    
    # tau
    v <- rinvgamma(1, 1, 1+1/tau^2)
    tau <- sqrt( rinvgamma(1, (p+1)/2, 0.5*sum(beta^2/(sig^2*lam[ID]*delta))+1/v) )
  }
  return(list(beta=beta.pos, sig=sig.pos))
}



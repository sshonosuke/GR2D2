###-------------------------------------------------------------###
###               Application to abalone data                   ###
###-------------------------------------------------------------###
rm(list=ls())

## R packages
library(ggplot2)
library(tidyr)
library(dplyr)
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


# evaluation points for covariates
m <- 20
X_origin_sub <- X_origin[sample(1:n, m),]
X_sub <- matrix(NA, m, p)
for(g in 1:G){
  X_sub[,ID==g] <- eval.basis(X_origin_sub[,g], Basis[[g]])
}


## MCMC
set.seed(1)
mc <- 3000
bn <- 1000

Fit <- list()
Fit[[1]] <- R2D2(Y, X, a=1/2, b=1/2, mc=mc, bn=bn) 
Fit[[2]] <- HS(Y, X, a=1/2, b=1/2, mc=mc, bn=bn) 
Fit[[3]] <- GIGG(Y, X, ID, a=1/p, b=1/2, mc=mc, bn=bn) 
Fit[[4]] <- GHS(Y, X, ID, mc=mc, bn=bn) 
Fit[[5]] <- R2D4_Laplace(Y, X, Group_ID=ID, phi_dist="Dir", mc=mc, bn=bn) 
Fit[[6]] <- R2D4_Laplace(Y, X, Group_ID=ID, phi_dist="LN", mc=mc, bn=bn) 


## posterior median
M <- length(Fit)
hBeta <- c()
for(k in 1:M){
  hBeta <- cbind(hBeta, apply(Fit[[k]]$beta, 2, median))
}


## Plot (coefficients)
df <- as.data.frame( log(1+abs(hBeta[,3:6])) )
names(df) <- c("GIGG", "GHS", "gR2D2-D", "gR2D2-L")
df$Index <- 1:nrow(df) 
df_long <- pivot_longer(df, cols=-Index, names_to="Pattern", values_to="Estimate")
vline_positions <- c(10.5, 20.5, 30.5, 40.5, 50.5, 60.5)


# Figure 2
pdf("Figure2.pdf", height=7, width=10)
ggplot(df_long, aes(x = Index, y = Estimate, color = Pattern, group = Pattern)) +
  geom_line(size=0.2, linetype="solid", alpha=0.8) +  
  geom_point(size=1, alpha=0.7) + 
  geom_vline(xintercept = vline_positions, linetype = "dashed", color = "gray30", alpha = 0.7) + 
  labs(title="", x="Index of coefficients", y="log(1 + Estimate)", color="Priors") +
  theme_minimal(base_size = 12) +  
  theme(
    plot.title=element_text(hjust=0.5), legend.position="right",
    panel.grid.major=element_line(color="gray95", size=0.2), 
    panel.grid.minor=element_line(color="gray95", size=0.2), 
  )
dev.off()



## Plot (regression function)
name <- c("Length", "Diameter", "Height", "Whole Weight", "Shucked Weight", "Viscera Weight", "Shell Weight")
plot_data <- list()
for(g in 1:7){
  x <- seq(min(X_origin[,g]), max(X_origin[,g]), length=100)
  x_eval <- eval.basis(x, Basis[[g]])
  y <- x_eval %*% hBeta[ID==g, 3:6]
  colnames(y) <- c("GIGG", "GHS", "gR2D2-D", "gR2D2-L")
  plot_data[[g]] <- data.frame(x=x, y, Variable=name[g], check.names=FALSE)
}

plot_data <- bind_rows(plot_data)
tidy_data <- pivot_longer(plot_data, cols=-c(x, Variable), names_to="Method", 
                          values_to="Value", names_transform=list(Method=as.character) )


# Figure 3
pdf("Figure3.pdf", height=7, width=10)
ggplot(tidy_data, aes(x=x, y=Value, color=Method)) +
  geom_line() +
  facet_wrap(~Variable, scales="free", ncol=4) +
  labs(x="", y="", color="Method") +
  theme_minimal() +
  theme(strip.text=element_text(size=12, face="bold"), legend.position="top")
dev.off()


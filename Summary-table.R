###-------------------------------------------------------------###
###                 Code for Tables 1 and 2                     ###
###-------------------------------------------------------------###
## The following files (output of 'Monte-Carlo-Simulation.R') are required:
# 'Sim-p50-result1.RData'
# 'Sim-p50-result2.RData'
# 'Sim-p50-result3.RData'
# 'Sim-p50-result4.RData'
# 'Sim-p50-result5.RData'
# 'Sim-p200-result1.RData'
# 'Sim-p200-result2.RData'
# 'Sim-p200-result3.RData'
# 'Sim-p200-result4.RData'
# 'Sim-p200-result5.RData'
# 'Sim-p1000-result1.RData'
# 'Sim-p1000-result2.RData'
# 'Sim-p1000-result3.RData'
# 'Sim-p1000-result4.RData'
# 'Sim-p1000-result5.RData'


## Main part
Method <- c("R2D2", "HS", "GIGG", "GHS", "gR2D2-D", "gR2D2-LN")
M <- length(Method)

S <- 5     # total number of scenarios 
n_vec <- rep(c(250, 150, 200), rep(2*M, 3))
p_vec <- rep(c(50, 200, 1000), rep(2*M, 3))
SNR_vec <- rep(c(rep(0.5, M), rep(0.7, M)), 3)

Tab1 <- Tab2 <- matrix(NA, 6*M, 2*S)
dimnames(Tab1)[[1]] <- dimnames(Tab2)[[1]] <- paste0("(n=", n_vec, ", p=", p_vec, ", SNR=", SNR_vec, ") ",  rep(Method, 6)) 
dimnames(Tab1)[[2]] <- paste0("Scenario ", c(1:5, 1:5), " (", rep(c("Null", "Non-null"), rep(5,2)), ")")
dimnames(Tab2)[[2]] <- paste0("Scenario ", c(1:5, 1:5), " (", rep(c("CP", "AL"), rep(5,2)), ")")

for(p in c(50, 200, 1000)){
  for(s in 1:5){
    load(paste0("Sim-p", p, "-result", s, ".RData"))
    # SNR=0.5
    Tab1[p_vec==p & SNR_vec==0.5, s] <- apply(na.omit(MSE_noise[,,1]), 2, mean)
    Tab1[p_vec==p & SNR_vec==0.5, S+s] <- apply(na.omit(MSE_signal[,,1]), 2, mean)
    Tab2[p_vec==p & SNR_vec==0.5, s] <- 100*apply(na.omit(CP_all[,,1]), 2, mean)
    Tab2[p_vec==p & SNR_vec==0.5, S+s] <- apply(na.omit(AL_all[,,1]), 2, mean)
    # SNR=0.7
    Tab1[p_vec==p & SNR_vec==0.7, s] <- apply(na.omit(MSE_noise[,,2]), 2, mean)
    Tab1[p_vec==p & SNR_vec==0.7, S+s] <- apply(na.omit(MSE_signal[,,2]), 2, mean)
    Tab2[p_vec==p & SNR_vec==0.7, s] <- 100*apply(na.omit(CP_all[,,2]), 2, mean)
    Tab2[p_vec==p & SNR_vec==0.7, S+s] <- apply(na.omit(AL_all[,,2]), 2, mean)
  }
}


# Output 
write.csv(round(Tab1, 2), file="Table1.csv")
write.csv(round(Tab2, 2), file="Table2.csv")

# setwd("C:/Users/baipl/Dropbox (UFL)/Cohort change points detection/code")
rm(list = ls())
library("sparsevar")
library("VARDetect")
library("MTS")
source("../../simulation_script.R")
source("../../auxiliary_functions.R")

#####################################################################################
### model basic parameters setting test
T <- 200
p <- 20
M <- 20
tau_true <- c(floor(T/3), floor(2*T/3))
tau_sub <- rep(0, M)
r_T <- 10
subjects <- vector('list', M)
set.seed(1)
for(i in 1:M){
    if(i <= 5){
        tau_sub[i] <- tau_true[1] + sample(seq(-r_T, r_T), 1)
    }else{
        tau_sub[i] <- tau_true[2] + sample(seq(-r_T, r_T), 1)
    }
}


### transition matrix generating
support_mat <- matrix(0, p, p)
for(k in 1:(p-1)){
    support_mat[k, k+1] <- 1
}
magnitudes <- sample(seq(0.55, 0.85, 0.1), M, replace = TRUE)
niter <- 50
cp_result <- vector('list', niter)
individual_result <- vector('list', niter)
cluster_result <- vector('list', niter)
for(epoch in 1:niter){
    ### generate subject time series
    for(i in 1:M){
        brks <- c(tau_sub[i], T+1)
        entry <- magnitudes[i]
        try <- simu_var(method = "sparse", 
                        nob = T, 
                        k = p,
                        lags = 1,
                        seed = epoch*i,
                        brk = brks,
                        sp_pattern = "off-diagonal", 
                        signals = c(-entry, entry))
        subjects[[i]] <- as.matrix(try$series)
    }
    
    ### apply the iterative clustering algorithm (ICA)
    thres <- 2*log(p)*log(T)
    fit <- iterative_selection(subjects, lambdas = c(0.1, 0.1), gamma = thres)
    cp_result[[epoch]] <- fit$final_cps
    individual_result[[epoch]] <- fit$each_cp
    cluster_result[[epoch]] <- fit$clusters
}

save(cp_result, file = "final_est_cp.RData")
save(individual_result, file = "individual_est_cps.RData")
save(cluster_result, file = "final_cluster.RData")


##################################################################################
### summary
niter <- 50
cp1 = cp2 <- c()
for(i in 1:niter){
    ncp <- length(cp_result[[i]])
    for(j in 1:ncp){
        if(cp_result[[i]][j] <= tau_true[1] + 15 && cp_result[[i]][j] >= tau_true[1] - 15){
            cp1 <- c(cp1, cp_result[[i]][j])
        }else if(cp_result[[i]][j] <= tau_true[2] + 15 && cp_result[[i]][j] >= tau_true[2] - 15){
            cp2 <- c(cp2, cp_result[[i]][j])
        }
    }
}
mean(cp1/T); sd(cp1/T)
mean(cp2/T); sd(cp2/T)

# calculate Hausdorff distance for the 
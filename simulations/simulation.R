rm(list = ls())
library("sparsevar")
library("VARDetect")
library("MTS")
source("simulation_script.R")
#####################################################################################
### model basic parameters setting test
T <- 200
p <- 20
M <- 10
tau_true <- floor(T/2)
tau_sub <- rep(0, M)
r_T <- 10
subjects <- vector('list', M)
set.seed(1)
for(i in 1:M){
    tau_sub[i] <- tau_true + sample(seq(-r_T, r_T), 1)
}

### transition matrix generating
support_mat <- matrix(0, p, p)
for(k in 1:(p-1)){
    support_mat[k, k+1] <- 1
}
magnitudes <- sample(seq(0.55, 0.85, 0.1), M, replace = TRUE)

niter <- 50
cp_result = runtimes <- rep(0, niter)
est_mats <- vector('list', niter)
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
    
    ### obtain the detected cps
    fit <- single.cp.detect(subjects, lags = c(1, 1), skip = 10)
    cp_result[epoch] <- fit$est.cp
    est_mats[[epoch]] <- cbind(fit$left_mats, fit$right_mats)
    runtimes[epoch] <- fit$runtime
    print(epoch)
}
save(cp_result, file = 'estimated_cp.RData')
save(est_mats, file = 'estimated_mats.RData')
save(runtimes, file = 'runtimes.RData')
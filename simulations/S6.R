# setwd("C:/Users/baipl/Dropbox (UFL)/Cohort change points detection/code")
rm(list = ls())
library("sparsevar")
library("VARDetect")
library("MTS")
library("mvtnorm")
source("../simulation_script.R")
#####################################################################################
### model basic parameters setting test
T <- 200
p <- 20
M <- 10
q.t <- 1
tau_true <- c(floor(T/2))
m0 <- length(tau_true)
m <- m0 + 1
tau_sub <- rep(0, M)
r_T <- 20
subjects <- vector('list', M)
set.seed(1)
for(i in 1:M){
    tau_sub[i] <- tau_true + sample(seq(-r_T, r_T), 1)
}

niter <- 50
cp_result <- rep(0, niter)
est_mats <- vector('list', niter)
runtimes <- rep(0, niter)
for(epoch in 1:niter){
    ### create random structure transition matrices and noise parameters
    e.sigma <- 0.1*diag(p)
    phi.full <- matrix(0, p, p*q.t*m)
    aa <- 0.8
    set.seed(epoch)
    for(mm in 1:m){
        phi.full[,((mm-1)*q.t*p+1):((mm)*q.t*p)] <- 0
        for (j in 1:(p-1)){
            bool_1 <- sample(0:2, 1, prob = c(0.1, 0.8, 0.1))
            x_shift = sample(0:4, 1)
            if (bool_1 > 0 && (j + x_shift[1:bool_1] <= p) ){
                phi.full[j,((mm-1)*q.t*p+j +  x_shift[1:bool_1])] <- -aa
            }
        }
        if(mm %% 2 == 0){
            phi.full[,((mm-1)*q.t*p+1):((mm)*q.t*p)] <- -phi.full[,((mm-1)*q.t*p+1):((mm)*q.t*p)]
        }
    }  
    
    ### generate subject time series
    for(i in 1:M){
        set.seed(epoch*i)
        brks <- c(tau_sub[i], T+1)
        try <- var.sim.break.tdist(nobs = T, 
                                   arlags = 1,
                                   phi = phi.full,
                                   sigma = e.sigma,
                                   brk = brks, 
                                   df = 5)
        subjects[[i]] <- as.matrix(try$series)
        # MTSplot(subjects[[i]])
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
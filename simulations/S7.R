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
M <- 20

tau_true <- floor(T/2)
m0 <- length(tau_true)
m <- m0 + 1

tau_sub <- rep(0, M)
r_T <- 10
subjects <- vector('list', M)
set.seed(1)
for(i in 1:M){
    tau_sub[i] <- tau_true + sample(seq(-r_T, r_T), 1)
}

### sample lags for each subjects
lag_sub <- rep(0, M)
set.seed(5)
for(i in 1:M){
    lag_sub[i] <- sample(c(1,2), 1)
}

### transition matrix signals generating
signal_sub <- vector('list', M)
set.seed(1)
for(i in 1:M){
    signal_sub[[i]] <- rep(0, lag_sub[i]*2)
    for(j in 1:(lag_sub[i]*2)){
        if(lag_sub[i] == 1){
            signal_sub[[i]][j] <- (-1)^j * sample(seq(0.5, 0.8, 0.1), 1)
        }else{
            signal_sub[[i]][j] <- (-1)^j * sample(seq(0.35, 0.55, 0.05), 1)
        }
    }
}

phi.full.subs <- vector('list', M)
for(i in 1:M){
    q.t <- lag_sub[i]
    temp <- matrix(0, p, p*q.t*m)
    
    for(j in 1:(p-1)){
        for(q in 1:(q.t*m)){
            aa <- signal_sub[[i]][q]
            temp[j, j+(q-1)*p+1] <- aa
        }
    }
    phi.full.subs[[i]] <- temp
}

e.sigma <- 0.1*diag(p)
niter <- 50
cp_result <- rep(0, niter)
est_mats <- vector('list', niter)
runtimes <- rep(0, niter)
for(epoch in 1:niter){
    ### generate subject time series
    for(i in 1:M){
        set.seed(epoch*i)
        q.t <- lag_sub[i]
        phi_full <- phi.full.subs[[i]]
        brks <- c(tau_sub[i], T+1)
        try <- var.sim.break.tdist(nobs = T, 
                                   arlags = 1:q.t, 
                                   phi = phi_full, 
                                   sigma = e.sigma,
                                   brk = brks,
                                   df = 10)
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
# setwd("~/Dropbox (UFL)/Cohort change points detection/code/simulations")
rm(list = ls())
library("sparsevar")
library("VARDetect")
library("cluster")
library("factoextra")
library("MTS")
source("../simulation_script.R")

### model basic parameters setting test
T <- 300
p <- 20
M <- 10
r_T <- 10
tau_true <- c(floor(T/3), floor(2*T/3))

tau_sub <- vector('list', M)
subjects <- vector('list', M)
set.seed(1)
for(i in 1:M){
    tau_sub[[i]] <- c(tau_true[1] + sample(seq(-r_T, r_T, 1), 1), 
                      tau_true[2] + sample(seq(-r_T, r_T, 1), 1))
}

#####################################################################
nepoch <- 50
skip <- 3
est_cp <- vector('list', nepoch)
est_mats <- vector('list', nepoch)
runtimes <- rep(0, nepoch)
for(epoch in 1:nepoch){
    ### transition matrix generating
    set.seed(epoch+9294)
    signals <- sample(seq(0.7, 0.85, 0.05), M, replace = TRUE)
    e.sigma <- 0.01*diag(p)
    
    ### subject time series generatings
    subjects <- vector('list', M)
    for(i in 1:M){
        brks <- c(tau_sub[[i]], T+1)
        entry <- signals[i]
        try <- simu_var(method = "sparse", 
                        nob = T, 
                        k = p,
                        sigma = e.sigma,
                        lags = 1,
                        seed = epoch*i,
                        brk = brks,
                        sp_pattern = "off-diagonal", 
                        signals = c(-entry, entry, -entry))
        subjects[[i]] <- as.matrix(try$series)
    }
    
    start_time <- proc.time()
    ### 1. detect change points for each subject
    cp_subjects <- c()
    for(i in 1:M){
        data <- subjects[[i]]
        fit <- tbss(data, lambda.2.cv = seq(0.05, 0.0005, length.out=10))
        cp_subjects <- c(cp_subjects, fit$cp)
        cat("\n")
        cat("Subject ", i, " ", "detected cps are: ", fit$cp, "\n")
    }
    
    ### 2. cluster the change points
    K <- get.optimal.K(cp_subjects)
    fit <- kmeans(cp_subjects, K)
    cp_cluster <- vector("list", K)
    for(k in 1:K){
        cp_cluster[[k]] <- cp_subjects[fit$cluster == k]
    }
    
    ### 3. exhaustive search on each cluster
    final.cp <- c()
    for(k in 1:K){
        cp_clust <- cp_cluster[[k]]
        lb <- min(cp_clust)
        ub <- max(cp_clust)
        if(ub - lb <= 2*skip+2){
            final.cp <- c(final.cp, floor((ub+lb)/2))
        }else{
            new_subjects <- vector('list', M)
            for(mm in 1:M){
                new_subjects[[mm]] <- subjects[[mm]][lb:ub, ]
            }
            exfit <- single.cp.detect(new_subjects, lags = c(1,1), skip = skip)
            final.cp <- c(final.cp, exfit$est.cp + skip + lb-1)
        }
    }
    cat("=======================================", "\n")
    cat(final.cp, "\n")
    cat("=======================================", "\n")
    end_time <- proc.time()
    
    ### 4. refit the lasso to obtain the estimated matrices for each subjects
    final.breaks <- c(1, final.cp, T)
    est_mat_subs <- vector('list', M)
    for(mm in 1:M){
        data <- subjects[[mm]]
        est.mat.subject <- NULL
        for(jj in 1:(length(final.breaks)-1)){
            s <- final.breaks[jj]
            e <- final.breaks[jj+1]
            data_segment <- data[s:e, ]
            fit_lasso <- fitVAR(data_segment, p = 1, nlambda = 5, nfolds = 5, 
                                threshold = TRUE)
            est.mat.subject <- cbind(est.mat.subject, fit_lasso$A[[1]])
        }
        est_mat_subs[[mm]] <- est.mat.subject
    }
    est_cp[[epoch]] <- final.cp
    est_mats[[epoch]] <- est_mat_subs
    runtimes[epoch] <- (end_time - start_time)[3]
    print(paste("finished epoch:", epoch, sep = " "))
    print("==============================================")
}

save(est_cp, file = "estimated_cps.RData")
save(est_mats, file = "estimated_mats.RData")
save(runtimes, file = "runningtimes.RData")

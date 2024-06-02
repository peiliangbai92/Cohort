# setwd("~/Dropbox (UFL)/Cohort change points detection/code/simulations")
rm(list = ls())
library("sparsevar")
library("VARDetect")
library("cluster")
library("factoextra")
library("mvtnorm")
library("MTS")
source("../simulation_script.R")

get.optimal.K.new <- function(cps){
    cps <- sort(cps)
    if(length(cps) >= 2){
        gap.temp <- sapply(2:length(cps), function(jjj) cps[jjj] - cps[jjj-1])
    }
    if(length(cps) > 5){
        if(length(unique(gap.temp)) > 1){
            print(fviz_nbclust(matrix(cps, length(cps), 1), kmeans, nstart = 25, method = "silhouette", 
                               k.max = min(5, length(cps)-1), nboot = 100) + 
                      labs(subtitle = "Gap statistic method"))
            cl <- fviz_nbclust(matrix(cps, length(cps), 1), kmeans, nstart = 25, method = "silhouette", 
                               k.max = min(5, length(cps)-1), nboot = 100) + 
                labs(subtitle = "Gap statistic method")
            cl.data <- cl$data
            gap <- cl.data$y
            cl.number <- which.max(diff(gap))+1
            gap.order <- order(gap, decreasing = TRUE)
        }
    }
    return(cl.number)
}

### model basic parameters setting test
T <- 300
p <- 20
M <- 10
r_T <- 10
q.t <- 1

tau_true <- c(floor(T/3), floor(2*T/3))
m0 <- length(tau_true)
m <- m0 + 1

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
    set.seed(epoch)
    
    ### create random structure transition matrices and noise parameters
    signals <- sample(seq(0.7, 0.85, 0.05), M, replace = TRUE)
    e.sigma <- 0.1*diag(p)
    phi.full <- matrix(0, p, p*q.t*m)
    aa <- 0.8
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
    for(i in 1:M){
        set.seed(epoch*i)
        brks <- c(tau_sub[[i]], T+1)
        try <- var.sim.break.tdist(nobs = T, 
                                   arlags = 1,
                                   phi = phi.full,
                                   sigma = e.sigma,
                                   brk = brks, 
                                   df = 4)
        subjects[[i]] <- as.matrix(try$series)
    }
    
    start_time <- proc.time()
    ### 1. detect change points for each subject
    cp_subjects <- c()
    for(i in 1:M){
        tryCatch({
            data <- subjects[[i]]
            fit <- tbss(data,
                        lambda.1.cv = seq(0.01, 0.0005, length.out=10),
                        lambda.2.cv = seq(0.025, 0.0005, length.out=10))
            cp_subjects <- c(cp_subjects, fit$cp)
            cat("\n")
            cat("Subject ", i, " ", "detected cps are: ", fit$cp, "\n")
        }, error = function(e){cat("ERROR:", conditionMessage(e), "\n")})
    }
    
    ### 2. cluster the change points
    K <- get.optimal.K.new(cp_subjects)
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
            tryCatch({
                new_subjects <- vector('list', M)
                for(mm in 1:M){
                    new_subjects[[mm]] <- subjects[[mm]][lb:ub, ]
                }
                exfit <- single.cp.detect(new_subjects, lags = c(1,1), skip = skip)
                final.cp <- c(final.cp, exfit$est.cp + skip + lb-1)
            }, error = function(e){cat("ERROR:", conditionMessage(e), "\n")})
        }
    }
    cat("Final cps are: ", final.cp, "\n")
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

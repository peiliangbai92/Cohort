### compare with SBS and BSS methods
setwd("~/Dropbox (UFL)/Cohort change points detection/code/comparison")

rm(list = ls())
library("sparsevar")
library("VARDetect")
library("factoextra")
library("cluster")
library("hdbinseg")
library("MTS")
source("../simulation_script.R")
source("../auxiliary_functions.R")

### model basic parameters setting test
T <- 400
p <- 20
M <- 20
tau_true <- c(floor(T/5), floor(T/2), floor(4*T/5))
tau_sub <- rep(0, M)
r_T <- 10
subjects <- vector('list', M)
set.seed(1)
for(i in 1:M){
    if(i <= 3){
        tau_sub[i] <- tau_true[1] + sample(seq(-r_T, r_T), 1)
    }else if(i <= 10){
        tau_sub[i] <- tau_true[2] + sample(seq(-r_T, r_T), 1)
    }else{
        tau_sub[i] <- tau_true[3] + sample(seq(-r_T, r_T), 1)
    }
}

### transition matrix generating
magnitudes <- sample(seq(0.55, 0.85, 0.1), M, replace = TRUE)
niter <- 50
cp_result <- vector('list', niter)
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
                        sp_pattern = "random", 
                        sp_density = c(0.05, 0.05),
                        signals = c(-entry, entry))
        print(brks)
        subjects[[i]] <- as.matrix(try$series)
    }
    
    ### apply the BSS method
    bss_cps <- vector('list', M)
    for(i in 1:M){
        fit_bss <- tryCatch(
            tbss(subjects[[i]], method = "sparse", 
                 lambda.1.cv = seq(0.001, 1, length.out=10), lambda.2.cv = seq(0.001, 1, length.out=10), 
                 block.size = 18),
            error = function(e){e}
        )
        cp <- fit_bss$cp
        if(length(cp) != 0){
            bss_cps[[i]] <- cp
        }
        print('\n')
        print(paste0("finished subject", i, sep = " "))
    }
    # deploy k-means for SBS results
    bss_cp_vec <- unlist(bss_cps)
    cl <- fviz_nbclust(matrix(bss_cp_vec, length(bss_cp_vec), 1), FUN = kmeans, nstart = 25, 
                       method = "gap_stat", k.max = 10, nboot = 100)
    cl.data <- cl$data
    gap <- cl.data$gap
    cur_gap <- gap[1]
    for(kk in 2:length(gap)){
        if(gap[kk] < cur_gap){
            cl.number <- kk-1
            break
        }else{
            cur_gap <- gap[kk]
        }
    }
    if(kk == length(gap)){
        cl.number <- length(gap)
    }
    bss_cl <- kmeans(bss_cp_vec, cl.number)
    
    bss_result <- cbind(bss_cp_vec, bss_cl$cluster)
    cluster_result[[epoch]] <- bss_result
}

save(cluster_result, file = "bss_cluster.RData")


### result analysis
load('bss_cluster.RData')
niter <- length(cluster_result)
nclusters <- rep(0, niter)
acc <- rep(NA, niter)
mean_grp1 = mean_grp2 = mean_grp3 <- rep(0, niter)
sd_grp1 = sd_grp2 = sd_grp3 <- rep(0, niter)
for(epoch in 1:niter){
    curr_result <- cluster_result[[epoch]]
    nclusters[epoch] <- max(curr_result[,2])
    
    # calculate the mis-classification rate
    group1 = group2 = group3 <- c()
    grp1_idx = grp2_idx = grp3_idx <- c()
    for(row in 1:dim(curr_result)[1]){
        # get estiamted cp information
        est_cp <- curr_result[row, 1]
        est_group <- curr_result[row, 2]
        
        # check if the current estimated cp belongs to the first true group
        if(est_cp < tau_true[1] + (tau_true[2] - tau_true[1]) / 2 && est_cp > tau_true[1] - (tau_true[1] - 1) / 2){
            group1 <- c(group1, est_cp)
            grp1_idx <- c(grp1_idx, est_group)
        }
        mean_grp1[epoch] <- mean(group1/T)
        sd_grp1[epoch] <- sd(group1/T)
        
        if(est_cp < tau_true[2] + (tau_true[3] - tau_true[2]) / 2 && est_cp > tau_true[2] - (tau_true[2] - tau_true[1]) / 2){
            group2 <- c(group2, est_cp)
            grp2_idx <- c(grp2_idx, est_group)
        }
        mean_grp2[epoch] <- mean(group2/T)
        sd_grp2[epoch] <- sd(group2/T)
        
        if(est_cp < tau_true[3] + (T - tau_true[3]) / 2 && est_cp > tau_true[3] - (tau_true[3] - tau_true[2]) / 2){
            group3 <- c(group3, est_cp)
            grp3_idx <- c(grp3_idx, est_group)
        }
        mean_grp3[epoch] <- mean(group3/T)
        sd_grp3[epoch] <- sd(group3/T)
    }
    
    # get the mis-classification error, first we need index for each class
    grp1 <- unique(grp1_idx)
    grp2 <- unique(grp2_idx)
    grp3 <- unique(grp3_idx)
    ttp = tfp = ttn = tfn <- 0
    
    # for group 1, we calcualte total accuracy
    for(row in 1:dim(curr_result)[1]){
        if(row <= 3){
            # group 1
            est_grp <- curr_result[row, 2]
            if(length(grp1) != 0){
                true_grp <- grp1
                if(est_grp == true_grp){
                    ttp <- ttp + 1
                }else{
                    tfn <- tfn + 1
                }
            }
        }else if(row > 3 && row <= 10){
            est_grp <- curr_result[row, 2]
            if(length(grp2) != 0){
                true_grp <- grp2
                if(est_grp == true_grp){
                    ttn <- ttn + 1
                }else{
                    tfp <- tfp + 1
                }
            }
        }else if(row > 10){
            est_grp <- curr_result[row, 2]
            if(length(grp3) != 0){
                true_grp <- grp3
                if(est_grp == true_grp){
                    ttn <- ttn + 1
                }else{
                    tfp <- tfp + 1
                }
            }
        }
    }
    
    # for group 2, we calcualte total accuracy
    for(row in 1:dim(curr_result)[1]){
        if(row <= 3){
            # group 1
            est_grp <- curr_result[row, 2]
            if(length(grp1) != 0){
                true_grp <- grp1
                if(est_grp == true_grp){
                    ttn <- ttn + 1
                }else{
                    tfp <- tfp + 1
                }
            }
            
        }else if(row > 3 && row <= 10){
            est_grp <- curr_result[row, 2]
            if(length(grp2) != 0){
                true_grp <- grp2
                if(est_grp == true_grp){
                    ttp <- ttp + 1
                }else{
                    tfn <- tfn + 1
                }
            }
            
        }else if(row > 10){
            est_grp <- curr_result[row, 2]
            if(length(grp3) != 0){
                true_grp <- grp3
                if(est_grp == true_grp){
                    ttn <- ttn + 1
                }else{
                    tfp <- tfp + 1
                }
            }
            
        }
    }
    
    # for group 3, we calcualte total accuracy
    for(row in 1:dim(curr_result)[1]){
        if(row <= 3){
            # group 1
            est_grp <- curr_result[row, 2]
            if(length(grp1) != 0){
                true_grp <- grp1
                if(est_grp == true_grp){
                    ttn <- ttn + 1
                }else{
                    tfp <- tfp + 1
                }
            }
            
        }else if(row > 3 && row <= 10){
            est_grp <- curr_result[row, 2]
            if(length(grp2) != 0){
                true_grp <- grp2
                if(est_grp == true_grp){
                    ttn <- ttn + 1
                }else{
                    tfp <- tfp + 1
                }
            }
            
        }else if(row > 10){
            est_grp <- curr_result[row, 2]
            if(length(grp3) != 0){
                true_grp <- grp3
                if(est_grp == true_grp){
                    ttp <- ttp + 1
                }else{
                    tfn <- tfn + 1
                }
            }
            
        }
    }
    acc[epoch] <- (ttp+ttn) / (ttp+ttn+tfp+tfn)
}

mean(nclusters)
mean(mean_grp1, na.rm = TRUE)
mean(sd_grp1, na.rm = TRUE)
mean(mean_grp2, na.rm = TRUE)
mean(sd_grp2, na.rm = TRUE)
mean(mean_grp3, na.rm = TRUE)
mean(sd_grp3, na.rm = TRUE)
mean(acc, na.rm = TRUE)


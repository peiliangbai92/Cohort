# setwd("C:/Users/baipl/Dropbox (UFL)/Cohort change points detection/code")

#####################################################################################
#' Generating non-stationary ARMA data.
#' 
#' @param nobs number of time points
#' @param arlags the true AR order
#' @param malags the true MA order
#' @param cnst the constant
#' @param phi parameter matrix of the AR model
#' @param theta parameter matrix of the MA model
#' @param sigma covariance matrix of the white noise
#' @return Matrices of time series data and white noise data 
var.sim.break.tdist <- function (nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL,
                                 theta = NULL, skip = 200, sigma, brk = nobs+1, df = 1 ) {
    if (!is.matrix(sigma)) 
        sigma = as.matrix(sigma)
    k <- nrow(sigma); m <- length(brk); nT <- nobs + skip
    
    #generate multivariate normal distributed data as the white noise data
    # at <- rmvnorm(nT, rep(0, k), sigma)
    at <- rmvt(nT, sigma = sigma, df = df)
    
    #generate the ARMA time series data
    nar <- length(arlags); p <- 0
    if (nar > 0) {
        arlags <- sort(arlags)
        p <- arlags[nar]
    }
    
    nma <- length(malags); q <- 0
    if (nma > 0) {
        malags <- sort(malags)
        q <- malags[nma]
    }
    
    ist = max(p, q) + 1
    zt = matrix(0, nT, k)
    if (length(cnst) == 0) 
        cnst = rep(0, k)
    if (m == 1){
        for (it in ist:nT) {
            tmp = matrix(at[it, ], 1, k)
            if (nma > 0) {
                for (j in 1:nma) {
                    jdx = (j - 1) * k
                    thej = theta[, (jdx + 1):(jdx + k)]
                    atm = matrix(at[it - malags[j], ], 1, k)
                    tmp = tmp - atm %*% t(thej)
                }
            }
            if (nar > 0) {
                for (i in 1:nar) {
                    idx = (i - 1) * k
                    phj = phi[, (idx + 1):(idx + k)]
                    ztm = matrix(zt[it - arlags[i], ], 1, k)
                    tmp = tmp + ztm %*% t(phj)
                }
            }
            zt[it, ] = cnst + tmp
        }
    }
    
    #if there are some break points
    if (m > 1){
        for (it in ist:(skip+brk[1]-1)) {
            tmp = matrix(at[it, ], 1, k)
            if (nma > 0) {
                for (j in 1:nma) {
                    jdx = (j - 1) * k
                    thej = theta[, (jdx + 1):(jdx + k)]
                    atm = matrix(at[it - malags[j], ], 1, k)
                    tmp = tmp - atm %*% t(thej)
                }
            }
            if (nar > 0) {
                for (i in 1:nar) {
                    idx = (i - 1) * k
                    phj = phi[, (idx + 1):(idx + k)]
                    ztm = matrix(zt[it - arlags[i], ], 1, k)
                    tmp = tmp + ztm %*% t(phj)
                }
            }
            zt[it, ] = cnst + tmp
        }
        for ( mm in 1:(m-1)){
            for (it in (skip+brk[mm]):(skip+brk[mm+1]-1) ) {
                tmp = matrix(at[it, ], 1, k)
                if (nma > 0) {
                    for (j in 1:nma) {
                        jdx = (j - 1) * k
                        thej = theta[, (jdx + 1):(jdx + k)]
                        atm = matrix(at[it - malags[j], ], 1, k)
                        tmp = tmp - atm %*% t(thej)
                    }
                }
                if (nar > 0) {
                    for (i in 1:nar) {
                        idx = (i - 1) * k
                        phj = phi[, ((mm)*p*k+idx + 1):((mm)*p*k+idx + k)]
                        ztm = matrix(zt[it - arlags[i], ], 1, k)
                        tmp = tmp + ztm %*% t(phj)
                    }
                }
                zt[it, ] = cnst + tmp
            }
        }
    }
    
    zt = zt[(1 + skip):nT, ]
    at = at[(1 + skip):nT, ]
    VARMAsim <- list(series = zt, noises = at)
}


### single change point selection procedure
single.cp.detect <- function(subjects, lags, skip = 5){
    n <- dim(subjects[[1]])[1]
    p <- dim(subjects[[1]])[2]
    M <- length(subjects)
    max.lag <- max(lags)
    
    # exhaustive search for change point
    obj_vals <- c()
    start <- proc.time()
    for(t in skip:(n-skip)){
        left_tmp = right_tmp <- 0
        for(i in 1:M){
            # separate data
            data <- subjects[[i]]
            left_seg <- data[1:t, ]
            right_seg <- data[(t+1):n, ]
            
            # fitting model parameters for both left and right sides
            fit_left <- fitVAR(left_seg, p = max.lag, nlambda = 5, 
                               nfolds = 5, threshold = TRUE)
            fit_right <- fitVAR(right_seg, p = max.lag, nlambda = 5, 
                                nfolds = 5, threshold = TRUE)
            
            # obtain the residuals and objective function
            left_tmp <- left_tmp + sum(fit_left$residuals^2)
            right_tmp <- right_tmp + sum(fit_right$residuals^2)
        }
        obj_vals <- c(obj_vals, left_tmp + right_tmp)
        # print(t)
    }
    end <- proc.time()
    elapsed_time <- (end - start)[3]
    idx <- which.min(obj_vals)
    cp <- idx + skip
    
    # refitting the left and right transition matrices
    left <- data[1:cp, ]; lag_left = lags[1]
    right <- data[(cp+1):n, ]; lag_right = lags[2]
    
    fit1 <- fitVAR(left, p = lag_left, threshold = TRUE)
    fit2 <- fitVAR(right, p = lag_right, threshold = TRUE)
    
    left_est <- fit1$A
    right_est <- fit2$A
    return(list(est.cp = cp, 
                left_mats = left_est,
                right_mats = right_est, 
                runtime = elapsed_time))
}


### multiple change points detection procedure
rolling.window.process <- function(subjects, lags, w, step.size, skip = 3){
    #'@param lags it should be a 2-element vector contains the lags for the left and right intervals, like c(1,1)
    #'@param w the window size
    #'@param step.size the rolling step size

    M <- length(subjects)
    nT <- dim(subjects[[1]])[1]; p <- dim(subjects[[1]])[2]
    s <- 1
    e <- s + w - 1
    candi <- c()
    while(e <= nT){
        window_subjects <- vector('list', M)
        for(i in 1:M){
            window_subjects[[i]] <- subjects[[i]][s:e, ]
        }
        ret <- single.cp.detect(window_subjects, lags = c(1, 1), skip = skip)
        cp <- ret$est.cp + s
        candi <- c(candi, cp)
        print("Finished:")
        cat(c(s, e, cp, '\n'))
        s <- s + step.size
        e <- e + step.size
    }
    return(unique(candi))
}

get.local.est.mat <- function(subjects, pts, an){
    M <- length(subjects)
    nT <- dim(subjects[[1]])[1]; p <- dim(subjects[[1]])[2]
    m <- length(pts)
    
    ### set up the boundary for the interval to compute the loss function of LIC
    bounds.1 <- vector('list', 2*m); bounds.2 <- vector('list', m)
    for(i in 1:m){
        if(pts[i] + an <= nT){
            bounds.1[[(2*i-1)]] <- c(pts[i] - an, pts[i] - 1)
            bounds.1[[(2*i)]] <- c(pts[i], pts[i] + an)
            bounds.2[[i]] <- c(pts[i] - an, pts[i] + an)
        }else{
            bounds.1[[(2*i-1)]] <- c(pts[i] - an, pts[i] - 1)
            bounds.1[[(2*i)]] <- c(pts[i], nT)
            bounds.2[[i]] <- c(pts[i] - an, nT)
        }
    }
    
    ### construct the vector for storing estimated matrices for each subject at each change point
    est.mat.subject <- vector('list', M)
    for(j in 1:M){
        est.mat.subject[[j]] <- vector('list', 2)
        est.mat.subject[[j]][[1]] <- vector('list', m)
        est.mat.subject[[j]][[2]] <- vector('list', m)
    }
    
    ### compute the loss function with the current change point
    L.n.1 <- c()
    for(mm in 1:(2*m)){
        tmp.sum <- 0
        for(j in 1:M){
            data <- subjects[[j]]
            data.temp <- data[(bounds.1[[mm]][1]):(bounds.1[[mm]][2]), ]
            try <- fitVAR(data.temp, nlambda = 5, nfolds = 5, threshold = TRUE)
            tmp.sum <- tmp.sum + sum(try$residuals^2)
            
            key <- ceiling(mm / 2)
            if(mm %% 2 == 1){
                est.mat.subject[[j]][[1]][[key]] <- try$A
            }else{
                est.mat.subject[[j]][[2]][[key+1]] <- try$A
            }
        }
        L.n.1 <- c(L.n.1, tmp.sum)
    }
    
    ### compute the loss function without the current change point
    L.n.2 <- c()
    for(mm in 1:m){
        tmp.sum <- 0
        for(j in 1:M){
            data <- subjects[[j]]
            data.temp <- data[(bounds.2[[m]][1]):(bounds.2[[m]][2]), ]
            try <- fitVAR(data.temp, nlambda = 5, nfolds = 5, threshold = TRUE)
            tmp.sum <- tmp.sum + sum(try$residuals^2)
        }
        L.n.2 <- c(L.n.2, tmp.sum)
    }
    
    return(list(L.n.1 = L.n.1, L.n.2 = L.n.2, est.mat.subject = est.mat.subject))
}

get.BIC <- function(subjects, q = 1, pts){
    M <- length(subjects)
    m <- length(pts)
    nT <- dim(subjects[[1]])[1]; p <- dim(subjects[[j]])[2]
    partitions <- c(1, pts, nT)
    
    ### construct a vectors to store the total BIC value
    bic.total <- c()
    
    ### update BIC scores for all subjects
    for(i in 1:(length(partitions)-1)){
        bic.temp <- 0
        s <- partitions[i]; e <- partitions[i+1]
        for(j in 1:M){
            data <- subjects[[j]]
            data.temp <- data[s:e, ]
            fit <- fitVAR(data.temp, q, nlambda = 5, nfolds = 5, thresholds = TRUE)
            sse <- sum(fit$residuals^2)
            df <- sum(fit$A[[1]] != 0)
            bic.temp <- bic.temp + ((e-s) * log(sse / (e-s)) + log(e-s) * df)
        }
        bic.total <- c(bic.total, bic.temp)
    }
    return(list(BIC = bic.total))
}

# test <- get.local.est.mat(subjects, candi, 10)

get.LIC <- function(subjects, q = 1, candidates, an){
    m <- length(candidates)
    M <- length(subjects)
    nT <- dim(subjects[[1]])[1]; p <- dim(subjects[[1]])[2]
    
    ### compute the local loss function LIC for each selected change point
    try <- get.local.est.mat(subjects, candidates, an)
    L.n.1 <- try$L.n.1
    L.n.2 <- try$L.n.2
    est_mat_subject <- try$est.mat.subject
    
    V <- rep(0, m)
    for(i in 1:m){
        V[i] <- (L.n.1[2*i-1] + L.n.1[2*i]) - L.n.2[i]
    }
    
    ### use BIC to determine the k-means
    BIC.diff <- 1
    BIC.old <- 1e8
    pts.select <- c()
    omega <- 0
    while(BIC.diff > 0 & length(unique(V)) > 1){
        pts.select.old <- pts.select
        omega.old <- omega
        
        ### use k-means to cluster vector v
        clus2 <- kmeans(V, 2)
        fit2 <- clus2$betweenss / clus2$totss
        if(fit2 < 0.20){
            omega <- max(V) + 1e-6
            pts.select <- c(pts.select)
            break
        }else{
            loc <- clus2$cluster
            if(clus2$centers[1] > clus2$centers[2]){
                ### case 1: we choose the first cluster as updated candidates
                omega <- min(V[which(loc == 1)]) - 1e-6
                loc.idx <- which(loc == 1)
            }else{
                ### case 2: we choose the second cluster as updated candidates
                omega <- min(V[which(loc == 2)]) - 1e-6
                loc.idx <- which(loc == 2)
            }
            pts.select <- sort(c(pts.select, candidates[loc.idx]))
            V[loc.idx] <- V[length(V)]
        }
        
        ### update the new candidate and their corresponding estimated transition matrices
        BIC.new <- sum(get.BIC(subjects, q, pts.select)$BIC)
        BIC.diff <- BIC.old - BIC.new
        BIC.old <- BIC.new
        if(BIC.diff <= 0){
            pts.select <- sort(pts.select.old)
            omega <- omega.old
            break
        }
    }
    
    ### select the final change points by LIC
    L.n.1.temp <- L.n.1; L.n.2.temp <- L.n.2
    L.n.plot <- rep(0, m+1)
    L.n.plot[1] <- sum(L.n.1) + m*omega
    mm <- 0
    lic <- 0
    add.temp <- 0
}

get.optimal.K <- function(cps){
    cps <- sort(cps)
    if(length(cps) >= 2){
        gap.temp <- sapply(2:length(cps), function(jjj) cps[jjj] - cps[jjj-1])
    }
    if(length(cps) > 5){
        if(length(unique(gap.temp)) > 1){
            print(fviz_nbclust(matrix(cps, length(cps), 1), kmeans, nstart = 25, method = "gap_stat",
                               k.max = min(10, length(cps)-1), nboot = 100) +
                      labs(subtitle = "Gap statistic method"))
            cl <- fviz_nbclust(matrix(cps, length(cps), 1), kmeans, nstart = 25, method = "gap_stat",
                               k.max = min(10, length(cps)-1), nboot = 100) +
                labs(subtitle = "Gap statistic method")
            cl.data <- cl$data
            gap <- cl.data$gap
            cl.number <- which.max(gap)
            gap.order <- order(gap, decreasing = TRUE)
        }
    }
    return(cl.number)
}

# LIC <- function(subjects, candidates, an){
#     M <- length(subjects)
#     nT <- dim(subjects[[1]])[1]
#     partitions <- c(1, candidates, nT)
#     lic = total_sse <- c()
#     est_mat_subject <- vector('list', M)
#     for(j in 1:M){
#         est_mat_subject[[j]] <- vector('list', length(partitions)-1)
#     }
#     for(i in 1:(length(partitions)-1)){
#         ### set up current partitioned data
#         s <- partitions[i] + an
#         e <- partitions[i+1] - an
#         segment_sse <- 0
#         for(j in 1:M){
#             data <- subjects[[j]]
#             curr_seg <- data[s:e, ]
#             fit <- fitVAR(curr_seg, nlambda = 5, nfolds = 5, threshold = TRUE)
#             est_mat_subject[[j]][[i]] <- fit$A[[1]]
#             segment_sse <- segment_sse + sum(fit$residuals^2)
#         }
#         total_sse <- c(total_sse, segment_sse)
#     }
#     for(i in 1:(length(partitions)-2)){
#         curr_cp <- partitions[i+1]
#         lic <- c(lic, total_sse[i] + total_sse[i+1])
#     }
#     return(list(LIC = lic, SSE = total_sse))
# }

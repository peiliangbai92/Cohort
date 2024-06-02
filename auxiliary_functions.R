### auxiliary functions
fista.sparse <- function(A, b, lambda, d, niter=100, backtracking = TRUE, phi.true=NULL){
    #'@description Estimating sparse matrix by using FISTA algorithm
    #'@param A design matrix
    #'@param b response vector
    #'@param lambda tuning parameter
    #'@param d dimension of design matrix
    #'@param niter maximum iterations for applying FISTA algorithm, default is 100
    #'@param backtracking a boolean variable, indicate if we use backtracking
    #'@param phi.true the true transition matirx, only available when simulation
    
    tnew = t <- 1
    x <- matrix(0, d, d)
    xnew <- x
    y <- x
    AtA <- t(A) %*% A
    Atb <- t(A) %*% b
    
    is.simu <- TRUE
    if(is.null(phi.true)) {
        is.simu <- FALSE
        phi.true <- diag(d)
    }
    
    obj.val = rel.err <- c()
    if(backtracking == TRUE){
        L <- norm(A, "2")^2 / 5
        eta <- 2
    }else{
        L <- norm(A, "2")^2
    }
    for(i in 1:niter){
        if(backtracking == TRUE){
            L.bar <- L
            flag <- FALSE
            while(flag == FALSE){
                prox <- prox.func(y, A, b, L.bar, lambda, AtA, Atb)
                if(f.func(prox, A, b) <= Q.func(prox, y, A, b, L.bar, AtA, Atb)){
                    flag <- TRUE
                }else{
                    L.bar <- L.bar * eta
                }
            }
            L <- L.bar
        }
        x <- xnew
        xnew <- prox
        t <- tnew
        tnew <- (1 + sqrt(1 + 4*t^2)) / 2
        y <- xnew + ((t - 1) / tnew) * (xnew - x)
        
        obj.val <- c(obj.val, f.func(xnew, A, b) + g.func(xnew, lambda))
        rel.err <- c(rel.err, norm(xnew - phi.true, "F") / norm(phi.true, "F"))
    }
    if(is.simu) {
        return(list(phi.hat = t(xnew), obj.vals = obj.val, rel.err = rel.err))
    }else {
        return(list(phi.hat = t(xnew), obj.vals = obj.val))
    }
}


shrinkage <- function(y, tau){
    #'@description Soft thresholding function for each vector
    #'@param y input vector
    #'@param tau threshold
    z <- matrix(0, nrow = nrow(y), ncol = ncol(y))
    for(i in 1:nrow(y)){
        for(j in 1:ncol(y)){
            z[i,j] <- sign(y[i,j]) * max(abs(y[i,j]) - tau, 0)
        }
    }
    return(z)
}


f.func <- function(x, A, b){
    #'@description Loss function
    #'@param x input matrix
    #'@param A transition matrix
    #'@param b response vector
    return((1/2) * norm(A %*% x - b, "F")^2)
}


gradf.func <- function(x, AtA, Atb){
    #'@description Gradient of loss function
    #'@param x input matrix
    #'@param AtA the result of A'A
    #'@param Atb the result of A'b
    return(AtA %*% x - Atb)
}


g.func <- function(x, lambda){
    #'@description sparse penalty
    #'@param x input vector
    #'@param lambda tuning parameter
    return(lambda*sum(x))
}


Q.func <- function(x, y, A, b, L, AtA, Atb){
    #'@description Q function, defined in FISTA algorithm
    #'@param x proximal result of transition matrix
    #'@param y transition matrix of previous iterative step
    #'@param A design matrix
    #'@param b response vector
    #'@param L learning rate
    #'@param AtA result of A'A
    #'@param Atb result of A'b
    return(f.func(y, A, b) + sum((x - y) * gradf.func(y, AtA, Atb)) + (1/2) * L * norm(x - y, "F")^2)
}


prox.func <- function(y, A, b, L, lambda, AtA, Atb){
    #'@description proximal function
    #'@param y the previous iteration of transition matrix
    #'@param A design matrix
    #'@param b response vector
    #'@param L learning rate
    #'@param lambda tuning parameter
    #'@param AtA result of A'A
    #'@param Atb result of A'b
    Y <- y - (1 / L) * gradf.func(y, AtA, Atb)
    return(shrinkage(Y, 2*lambda / L))
}


hausdorff_dist <- function(A, B) {
    n <- length(A)
    m <- length(B)
    dist_mat <- matrix(0, nrow=n, ncol=m)
    for(r in 1:n) {
        for(c in 1:m) {
            dist_mat[r,c] <- abs(A[r] - B[c])
        }
    }
    min_each_row <- apply(dist_mat, MARGIN = 1, FUN='min')
    return(max(min_each_row))
}

L <- function(data, lambda, scale=2) {
    # objective function
    n <- dim(data)[1]
    p <- dim(data)[2]
    
    if(n < scale*p+1) {
        print("not enough time points")
        break
    }
    sse <- c()
    for(t in (scale*p+1):(n-(scale*p+1))) {
        left_data <- data[1:t,]
        right_data <- data[(t+1):n, ]
        
        # partition the left and right data
        n_left <- dim(left_data)[1]
        n_right <- dim(right_data)[1]
        X_left <- left_data[1:(n_left-1), ]
        Y_left <- left_data[2:n_left, ]
        X_right <- right_data[1:(n_right-1), ]
        Y_right <- right_data[2:n_right, ]
        
        if(n_left <= n_right) {
            left_fit <- fista.sparse(X_left, Y_left, lambda[1], p, niter=20)
            right_fit <- fista.sparse(X_right, Y_right, lambda[2], p, niter=20)
        }else {
            left_fit <- fista.sparse(X_left, Y_left, lambda[2], p, niter=20)
            right_fit <- fista.sparse(X_right, Y_right, lambda[1], p, niter=20)
        }
        
        # estimated transition matrices
        left_mat <- left_fit$phi.hat
        right_mat <- right_fit$phi.hat
        
        # get the estimated transition matrices
        left_res <- data[2:t,] - data[1:(t-1),] %*% t(left_mat)
        right_res <- data[(t+2):n,] - data[(t+1):(n-1),] %*% t(right_mat)
        
        left_sse <- sum(left_res^2)
        right_sse <- sum(right_res^2)
        
        sse <- c(sse, left_sse + right_sse)
        # print(t)
    }
    return(sse)
}


### iterative clustering algorithm (ICA) to detect abnormal subjects
iterative_selection <- function(subjects, gamma=NULL, lambdas, scale=2, refit=TRUE) {
    #'@param subjects a list of data for all subjects
    #'@param threshold a threshold for clustering
    
    # check if the list is empty
    N <- length(subjects)
    if (N == 0) {
        print("Empty list!")
        break
    }
    
    # record the dimension of the subject data
    temp_ts <- subjects[[1]]
    n <- dim(temp_ts)[1]
    p <- dim(temp_ts)[2]
    
    if(is.null(gamma)) {
        gamma <- 3*log(p)*log(n)
    }
    
    ###### Initialization step: find CPs for each subject by using BSS
    cps <- rep(0, N)
    subjects_sse <- vector('list', N)
    clusters <- vector('list', N)
    for(i in 1:N) {
        ts <- subjects[[i]]
        sse <- L(ts, lambda=lambdas, scale=scale)
        
        subjects_sse[[i]] <- sse
        cps[i] <- which.min(sse)+scale*p+1
        clusters[[i]] <- i
        print(paste('estimated CP is:', cps[i], sep = " "))
        print("============================================")
    }
    
    remain <- clusters
    nclusters <- 1
    final_clusters <- vector('list', N)
    while(length(remain) > 0) {
        nn <- length(remain)
        if (nn == 1) {
            final_clusters[[nclusters]] <- remain[[1]]
            remain[[1]] <- NULL
            break
        }
        distants <- rep(0, nn-1)
        for(j in 2:nn) {
            grp_1 <- remain[[1]]
            grp_2 <- remain[[j]]
            sse1 = sse2 <- 0
            for(ix in grp_1) {
                sse1 <- sse1 + subjects_sse[[ix]]
            }
            for(ix in grp_2) {
                sse2 <- sse2 + subjects_sse[[ix]]
            }
            
            # get cp from group 1 and group 2 separately
            cp_from_grp1 <- which.min(sse1) + (scale*p+1)
            cp_from_grp2 <- which.min(sse2) + (scale*p+1)
            
            # estimate cp by combining group 1 and 2
            cp_from_combine <- which.min((sse1+sse2) / (length(grp_1)+length(grp_2))) + (scale*p+1)
            
            # calculate Hausdorff distance between these two sets
            A <- c(cp_from_grp1, cp_from_grp2)
            B <- c(cp_from_combine)
            distants[j-1] <- hausdorff_dist(A, B)
        }
        
        # combine the minimum distance to the first group
        dmin <- min(distants)
        min_idx <- which.min(distants)+1
        if(dmin < gamma) {
            min_grp <- remain[[min_idx]]
            remain[[min_idx]] <- NULL
            remain[[1]] <- c(remain[[1]], min_grp)
        }else {
            final_clusters[[nclusters]] <- remain[[1]]
            remain[[1]] <- NULL
            nclusters <- nclusters + 1
        }
    }
    results <- final_clusters[1:nclusters]
    print(results)
    
    ### refit the change points for each cluster
    final_cps <- rep(0, nclusters)
    for(i in 1:nclusters) {
        grp <- final_clusters[[i]]
        sse <- 0
        for(sub in grp) {
            sse <- sse + subjects_sse[[sub]]
        }
        sse <- sse / length(grp)
        final_cps[i] <- which.min(sse) + (scale*p+1)
    }
    return(list(final_cps = final_cps, each_cp = cps,
                clusters = results, sse_list = subjects_sse))
}

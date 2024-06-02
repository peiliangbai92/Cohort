setwd("~/Dropbox (UFL)/Cohort change points detection/code/simulations")
rm(list = ls())
library("sparsevar")
library("VARDetect")
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

### transition matrix generating
epoch <- 1
set.seed(1)
signals <- sample(seq(0.7, 0.85, 0.05), M, replace = TRUE)
e.sigma <- 0.01*diag(p)
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

#####################################################################################
### step 1: rolling window step:
# w <- 20
# s <- 1
# e <- s + w - 1
# step.size <- 10
# candi <- c()
# while(e <= T){
#     window_subjects <- vector('list', M)
#     for(i in 1:M){
#         window_subjects[[i]] <- subjects[[i]][s:e, ]
#     }
#     ret <- single.cp.detect(window_subjects, lags = c(1, 1), skip = 3)
#     cp <- ret$est.cp + s
#     candi <- c(candi, cp)
#     print("Finished:")
#     cat(c(s, e, cp, '\n'))
#     s <- s + step.size
#     e <- e + step.size
# }
candidates <- rolling.window.process(subjects, lags = c(1, 1), w = 60, step.size = 25, skip = 10)

### step 2: calculate LIC for each candidate
get_lic <- LIC(subjects, candi)
lic <- get_lic$LIC
init_sse <- get_lic$SSE

# k-means + BIC step
fit <- kmeans(lic, 2)
cluster1 <- partitions[fit$cluster == 1]
cluster2 <- partitions[fit$cluster == 2]
lic_c1 <- lic[fit$cluster == 1]
lic_c2 <- lic[fit$cluster == 2]

large_set = small_set <- NULL
if(mean(lic_c1) < mean(lic_c2)){
    large_set <- 1
    small_set <- 2
    thres <- (min(lic_c1) + max(lic_c2)) / 2
}else{
    large_set <- 2
    small_set <- 1
    thres <- (min(lic_c2) + max(lic_c1)) / 2
}

bic_old <- 1e3
bic_new <- 0
J <- partitions[fit$cluster == large_set]
while(bic_new - bic_old < 0){
    bic_new
}
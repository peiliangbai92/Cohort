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
p <- 80
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

load("estimated_cps.RData")
load("estimated_mats.RData")
load("runningtimes.RData")

niter <- 50
cp1 = cp2 <- c()
for(i in 1:niter){
    cps <- est_cp[[i]]
    for(j in 1:length(cps)){
        if(cps[j] <= 120 && cps[j] >= 80){
            cp1 <- c(cp1, cps[j])
        }else if(cps[j] <= 220 && cps[j] >= 180){
            cp2 <- c(cp2, cps[j])
        }
    }
}
cat("selection rate for 1st cp is:", length(cp1) / niter, "\n")
cat("selection rate for 2nd cp is:", length(cp2) / niter, "\n")

mean(cp1/T); sd(cp1/T)
mean(cp2/T); sd(cp2/T)
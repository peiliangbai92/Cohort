### Comparison method: factor-adjusted VAR segmentation (fvarseg::fvar.seg)
### + kmeans (gap-stat) on G4 setting.
###
### Note: fvar.seg expects the input matrix in (variable x time) layout,
### i.e. p rows and T columns -- this is the convention used by its own
### sim.data() helper. simu_var() from VARDetect returns (time x variable),
### so we transpose before calling fvar.seg.

rm(list = ls())
suppressPackageStartupMessages({
    library("sparsevar")
    library("VARDetect")
    library("factoextra")
    library("fvarseg")
    library("MTS")
})
source("../simulation_script.R")

## --- G4 setting -----------------------------------------------------------
T <- 400; p <- 20; M <- 20
tau_true <- c(floor(T/8), floor(T/2), floor(4*T/5))
group_sizes <- c(3, 7, 10)
r_T <- 10
niter <- 50
output_file <- "compare_fvarseg.RData"

set.seed(1)
tau_sub <- integer(M)
true_label <- rep(seq_along(group_sizes), group_sizes)
for (i in 1:M) tau_sub[i] <- tau_true[true_label[i]] + sample(seq(-r_T, r_T), 1)
magnitudes <- sample(seq(0.55, 0.85, 0.1), M, replace = TRUE)

## --- main loop -----------------------------------------------------------
all_cps   <- vector("list", niter)
all_clust <- vector("list", niter)
all_subj_id <- vector("list", niter)
runtimes <- numeric(niter)

for (epoch in 1:niter) {
    subjects <- vector("list", M)
    for (i in 1:M) {
        brks <- c(tau_sub[i], T + 1)
        sim <- simu_var(method = "sparse", nob = T, k = p, lags = 1,
                        seed = epoch * i, brk = brks,
                        sp_pattern = "random", sp_density = c(0.05, 0.05),
                        signals = c(-magnitudes[i], magnitudes[i]))
        subjects[[i]] <- as.matrix(sim$series)
    }

    t0 <- proc.time()
    fv_cps <- vector("list", M)
    for (i in 1:M) {
        fit <- tryCatch(
            ## fvar.seg expects p x T; simu_var returns T x p.
            ## CV path slows the call by ~15x without changing the
            ## detection story for T=400, so we rely on the default
            ## fvar.seg threshold.
            fvar.seg(t(subjects[[i]])),
            error = function(e) list(common.out = list(est.cp = integer(0)),
                                     idio.out   = list(est.cp = integer(0))))
        cp <- sort(unique(c(fit$common.out$est.cp, fit$idio.out$est.cp)))
        fv_cps[[i]] <- cp
    }

    cp_vec <- numeric(0); sub_id_vec <- integer(0)
    for (i in 1:M) {
        if (length(fv_cps[[i]]) > 0) {
            cp_vec     <- c(cp_vec, fv_cps[[i]])
            sub_id_vec <- c(sub_id_vec, rep(i, length(fv_cps[[i]])))
        }
    }

    cl_labels <- integer(length(cp_vec))
    if (length(cp_vec) > 1) {
        cl <- tryCatch(
            fviz_nbclust(matrix(cp_vec, ncol = 1), FUN = kmeans,
                         nstart = 25, method = "gap_stat",
                         k.max = min(6, length(cp_vec) - 1), nboot = 25),
            error = function(e) NULL)
        if (!is.null(cl)) {
            gap <- cl$data$gap
            cl.number <- which.max(gap[-1] - gap[-length(gap)] < 0)
            if (length(cl.number) == 0 || cl.number < 1) cl.number <- 1
            cl_labels <- kmeans(cp_vec, cl.number)$cluster
        }
    }

    all_cps[[epoch]]   <- cp_vec
    all_clust[[epoch]] <- cl_labels
    all_subj_id[[epoch]] <- sub_id_vec
    runtimes[epoch]    <- (proc.time() - t0)[3]
    cat(sprintf("epoch %d done, |cps|=%d, K_hat=%d, t=%.1fs\n",
                epoch, length(cp_vec), max(c(cl_labels, 0)), runtimes[epoch]))
    save(all_cps, all_clust, all_subj_id, runtimes,
         tau_true, group_sizes, T, M, r_T,
         file = output_file)
}
cat("Saved to", output_file, "\n")

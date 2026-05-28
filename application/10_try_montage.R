### Re-run the whole application pipeline with an alternative 21-channel
### montage and compare Hamming patterns against the main (10-20 + mastoids)
### run. Outputs go to application/cache_B/ and application/output_B/.
###
### The goal is to find whether a different 21-channel subset recovers the
### "within-cluster Hamming < between-cluster Hamming" pattern reported in
### the paper but only partially reproduced under the 10-20 + mastoids set.

suppressPackageStartupMessages({
    library(edfReader); library(signal); library(igraph); library(pheatmap)
})

source("auxiliary_functions.R")
source("application/01_preprocess_bdf.R")        # imports robust_detrend, alpha_bandpass
source("application/02_build_intervals.R")       # imports window_subject
source("application/03_cohort_cp.R")             # imports fit_one_interval (uses sources above)
source("application/04_ica_clustering.R")        # imports cluster_from_sse
source("application/05_networks.R")              # imports .fit_segment, .threshold etc.
source("application/06_hamming.R")               # imports compute_hamming_tables

## ---- Alternative montage (Subset B): 19 standard 10-20 + POz + Oz, no mastoids
TARGET_CHANNELS_B <- c(
    "Fp1","Fp2","F7","F3","Fz","F4","F8",
    "T7", "C3","Cz","C4","T8",
    "P7", "P3","Pz","P4","P8",
    "POz","Oz",
    "O1", "O2"
)
stopifnot(length(TARGET_CHANNELS_B) == 21)

CACHE_B <- "application/cache_B"
OUT_B   <- "application/output_B"
dir.create(CACHE_B, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_B,   showWarnings = FALSE, recursive = TRUE)


## ---- Step 1: preprocess each BDF into a 21-channel matrix using Subset B
preprocess_B <- function(force = FALSE) {
    out <- vector("list", 22)
    for (sid in 1:22) {
        cache <- file.path(CACHE_B, sprintf("preproc_S%02d.rds", sid))
        if (!force && file.exists(cache)) {
            out[[sid]] <- readRDS(cache)
            cat(sprintf("[%2d/22] cached -> %s\n", sid, cache))
            next
        }
        bdf <- sprintf("raw_data/EEG_Cat_Study4_Resting_S%d.bdf", sid)
        t0  <- proc.time()
        X   <- preprocess_one_bdf(bdf, channels = TARGET_CHANNELS_B)
        saveRDS(X, cache); out[[sid]] <- X
        cat(sprintf("[%2d/22] %s  (%.1fs, %d x %d)\n",
                    sid, basename(bdf), (proc.time() - t0)[3], nrow(X), ncol(X)))
    }
    invisible(out)
}

## ---- Step 2: window into 3 intervals around the empirical CPs
build_intervals_B <- function(preproc) {
    intervals <- vector("list", length(EMPIRICAL_CPS))
    for (k in seq_along(EMPIRICAL_CPS)) {
        intervals[[k]] <- lapply(preproc, window_subject, cp = EMPIRICAL_CPS[k])
    }
    saveRDS(intervals, file.path(CACHE_B, "intervals.rds"))
    intervals
}

## ---- Step 3: cohort CP detection (slow step, ~40 min)
run_cohort_B <- function(intervals, force = FALSE) {
    sse_file <- file.path(CACHE_B, "sse_lists.rds")
    if (!force && file.exists(sse_file)) {
        cat("[cohort_B] cached -> ", sse_file, "\n"); return(readRDS(sse_file))
    }
    fits <- vector("list", length(intervals))
    for (k in seq_along(intervals)) {
        cat(sprintf("\n=== Subset B  Interval %d ===\n", k))
        fits[[k]] <- fit_one_interval(intervals[[k]])
        cat(sprintf("  cohort cp = %d   per-subject mean = %.1f   sd = %.1f\n",
                    fits[[k]]$cohort_cp,
                    mean(fits[[k]]$per_sub_cp), sd(fits[[k]]$per_sub_cp)))
    }
    saveRDS(fits, sse_file); invisible(fits)
}


## ---- Driver
if (sys.nframe() == 0L) {
    cat("\n>>> Using alternative montage (Subset B): no mastoids, +POz/Oz\n\n")
    preproc   <- preprocess_B()
    intervals <- build_intervals_B(preproc)
    fits      <- run_cohort_B(intervals)

    p <- fits[[1]]$p; n <- fits[[1]]$n
    cat("\n--- ICA clustering (gammas tuned same as main run) ---\n")
    mild <- lapply(fits, function(f)
        cluster_from_sse(f$sse_list, p, n, gamma = 200))
    cons <- lapply(fits, function(f)
        cluster_from_sse(f$sse_list, p, n, gamma = 375))
    for (k in 1:3) cat(sprintf("  CP%d  mild #clust = %d   cons #clust = %d\n",
                               k, length(mild[[k]]$clusters), length(cons[[k]]$clusters)))

    cat("\n--- Per-subject sparse VAR(1) networks ---\n")
    per_subj <- vector("list", 3)
    for (k in 1:3) {
        cps <- fits[[k]]$per_sub_cp
        L_list <- vector("list", length(intervals[[k]]))
        R_list <- vector("list", length(intervals[[k]]))
        for (i in seq_along(intervals[[k]])) {
            X <- intervals[[k]][[i]]; cp <- cps[i]
            n_L <- cp; n_R <- nrow(X) - cp
            lam_L <- if (n_L <= n_R) LAMBDAS[1] else LAMBDAS[2]
            lam_R <- if (n_R <  n_L) LAMBDAS[1] else LAMBDAS[2]
            L_list[[i]] <- .fit_segment(X[1:cp, , drop = FALSE], lam_L)
            R_list[[i]] <- .fit_segment(X[(cp + 1):nrow(X), , drop = FALSE], lam_R)
        }
        per_subj[[k]] <- list(left = L_list, right = R_list)
    }
    saveRDS(per_subj, file.path(CACHE_B, "per_subject_A.rds"))

    cat("\n--- Hamming distances under Subset B ---\n")
    ica_B <- list(mild = mild, conservative = cons)
    H <- compute_hamming_tables(per_subj, ica_B)
    for (setting in c("mild", "conservative")) for (k in 1:3) {
        M <- round(H[[setting]][[k]]$H_avg, 3)
        cat(sprintf("\nSubset B  %s  CP%d  (thr=%.3f):\n",
                    setting, k, H[[setting]][[k]]$thr))
        print(M)
        write.csv(M, file.path(OUT_B, sprintf("hamming_%s_cp%d.csv", setting, k)))
    }
    saveRDS(H, file.path(CACHE_B, "hamming.rds"))

    cat("\n--- Comparison: within vs max between Hamming, per cluster ---\n")
    for (setting in c("mild", "conservative")) for (k in 1:3) {
        M <- H[[setting]][[k]]$H_avg
        diag_vals <- diag(M); offd <- M; diag(offd) <- NA
        max_offd  <- apply(offd, 1, max, na.rm = TRUE)
        viol <- sum(diag_vals > max_offd)
        cat(sprintf("  %-12s CP%d:  diag = [%s]  max off-diag per row = [%s]  violators = %d/%d\n",
                    setting, k,
                    paste(round(diag_vals, 3), collapse = ", "),
                    paste(round(max_offd,  3), collapse = ", "),
                    viol, length(diag_vals)))
    }
}

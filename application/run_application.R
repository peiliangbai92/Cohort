### Re-run the application pipeline using the exact 21-channel montage from
### the original Figure 3 of the paper:
###
###   row 1 (prefrontal):  FP1   FPZ   FP2
###   row 2 (frontal):     F7  F3  FZ  F4  F8
###   row 3 (central):     T7  C3  CZ  C4  T8
###   row 4 (parietal):    P7  P3  PZ  P4  P8
###   row 5 (occipital):   O1   OZ   O2
###
### This is the 19 standard 10-20 channels plus the midline FPZ and OZ, with
### NO mastoids (the previous run used 19 + M1/M2 which gave a worse Hamming
### within < between pattern).
###
### All intermediates land in application/cache/, hamming tables in
### application/output/, and the publication PDFs used by the manuscript in
### tex/Cohort_cp_detection_manuscript/images/.

suppressPackageStartupMessages({
    library(edfReader); library(signal); library(igraph); library(pheatmap)
})

source("auxiliary_functions.R")
source("application/01_preprocess_bdf.R")
source("application/02_build_intervals.R")
source("application/03_cohort_cp.R")
source("application/04_ica_clustering.R")
source("application/05_networks.R")
source("application/06_hamming.R")

## Replace auxiliary_functions::fista.sparse with a fixed-Lipschitz variant.
## The backtracking branch in the original FISTA doubles L every time the Q
## condition fails; on poorly-conditioned EEG segments L explodes and step
## size collapses, which is why some subjects took >4000 s in the previous
## run. The non-backtracking variant uses L = ||A||^2 (the true Lipschitz
## bound), avoiding the pathology while keeping the same convergence
## guarantees. Each L() call drops from ~30 s (nominal) / ~4000 s (worst case)
## to a flat ~17 s.
fista.sparse <- function(A, b, lambda, d, niter = 20,
                         backtracking = TRUE, phi.true = NULL) {
    tnew <- t <- 1
    x <- xnew <- y <- matrix(0, d, d)
    AtA <- crossprod(A); Atb <- crossprod(A, b)
    L_lip <- norm(A, "2")^2
    for (i in seq_len(niter)) {
        prox <- prox.func(y, A, b, L_lip, lambda, AtA, Atb)
        x <- xnew; xnew <- prox
        t <- tnew; tnew <- (1 + sqrt(1 + 4 * t^2)) / 2
        y <- xnew + ((t - 1) / tnew) * (xnew - x)
    }
    list(phi.hat = t(xnew))
}

TARGET_CHANNELS_V2 <- c(
    "Fp1", "Fpz", "Fp2",
    "F7",  "F3",  "Fz",  "F4", "F8",
    "T7",  "C3",  "Cz",  "C4", "T8",
    "P7",  "P3",  "Pz",  "P4", "P8",
    "O1",  "Oz",  "O2"
)
stopifnot(length(TARGET_CHANNELS_V2) == 21)

CACHE_V2 <- "application/cache"
OUT_V2   <- "application/output"
IMG_V2   <- "tex/Cohort_cp_detection_manuscript/images"
for (d in c(CACHE_V2, OUT_V2, IMG_V2)) dir.create(d, showWarnings = FALSE, recursive = TRUE)


## ============= preprocessing =============================================
preproc_v2 <- function(force = FALSE) {
    out <- vector("list", 22)
    for (sid in 1:22) {
        cache <- file.path(CACHE_V2, sprintf("preproc_S%02d.rds", sid))
        if (!force && file.exists(cache)) {
            out[[sid]] <- readRDS(cache)
            cat(sprintf("[%2d/22] cached -> %s\n", sid, cache)); next
        }
        bdf <- sprintf("raw_data/EEG_Cat_Study4_Resting_S%d.bdf", sid)
        t0  <- proc.time()
        X   <- preprocess_one_bdf(bdf, channels = TARGET_CHANNELS_V2)
        saveRDS(X, cache); out[[sid]] <- X
        cat(sprintf("[%2d/22] %s  (%.1fs)\n",
                    sid, basename(bdf), (proc.time() - t0)[3]))
    }
    invisible(out)
}


intervals_v2 <- function(preproc, force = FALSE) {
    f <- file.path(CACHE_V2, "intervals.rds")
    if (!force && file.exists(f)) { cat("[intervals_v2] cached\n"); return(readRDS(f)) }
    intervals <- vector("list", length(EMPIRICAL_CPS))
    for (k in seq_along(EMPIRICAL_CPS))
        intervals[[k]] <- lapply(preproc, window_subject, cp = EMPIRICAL_CPS[k])
    saveRDS(intervals, f); invisible(intervals)
}


cohort_v2 <- function(intervals, force = FALSE) {
    f <- file.path(CACHE_V2, "sse_lists.rds")
    if (!force && file.exists(f)) { cat("[cohort_v2] cached\n"); return(readRDS(f)) }
    fits <- vector("list", length(intervals))
    for (k in seq_along(intervals)) {
        cat(sprintf("\n=== v2  Interval %d ===\n", k))
        fits[[k]] <- fit_one_interval(intervals[[k]])
        cat(sprintf("  cohort cp = %d  mean = %.1f  sd = %.1f\n",
                    fits[[k]]$cohort_cp,
                    mean(fits[[k]]$per_sub_cp), sd(fits[[k]]$per_sub_cp)))
    }
    saveRDS(fits, f); invisible(fits)
}


networks_v2 <- function(intervals, fits, force = FALSE) {
    f <- file.path(CACHE_V2, "per_subject_A.rds")
    if (!force && file.exists(f)) { cat("[networks_v2] cached\n"); return(readRDS(f)) }
    out <- vector("list", 3)
    for (k in 1:3) {
        cps <- fits[[k]]$per_sub_cp; M <- length(intervals[[k]])
        L_list <- R_list <- vector("list", M)
        for (i in 1:M) {
            X <- intervals[[k]][[i]]; cp <- cps[i]
            n_L <- cp; n_R <- nrow(X) - cp
            lam_L <- if (n_L <= n_R) LAMBDAS[1] else LAMBDAS[2]
            lam_R <- if (n_R <  n_L) LAMBDAS[1] else LAMBDAS[2]
            L_list[[i]] <- .fit_segment(X[1:cp, , drop = FALSE], lam_L)
            R_list[[i]] <- .fit_segment(X[(cp + 1):nrow(X), , drop = FALSE], lam_R)
        }
        out[[k]] <- list(left = L_list, right = R_list)
    }
    saveRDS(out, f); invisible(out)
}


## ============= plotting helpers ===========================================

## Manual hierarchical layout matching the paper's Fig 3 row structure:
##   row 1 (y=5):  FP1   FPZ   FP2
##   row 2 (y=4):  F7  F3  FZ  F4  F8
##   row 3 (y=3):  T7  C3  CZ  C4  T8
##   row 4 (y=2):  P7  P3  PZ  P4  P8
##   row 5 (y=1):  O1   OZ   O2
LAYOUT_V2 <- (function() {
    rows <- list(
        c("Fp1", "Fpz", "Fp2"),
        c("F7", "F3", "Fz", "F4", "F8"),
        c("T7", "C3", "Cz", "C4", "T8"),
        c("P7", "P3", "Pz", "P4", "P8"),
        c("O1", "Oz", "O2")
    )
    coords <- matrix(0, nrow = 21, ncol = 2,
                     dimnames = list(TARGET_CHANNELS_V2, c("x", "y")))
    for (r in seq_along(rows)) {
        ch <- rows[[r]]; nC <- length(ch)
        x  <- seq(-(nC - 1) / 2, (nC - 1) / 2, length.out = nC)
        coords[ch, "x"] <- x
        coords[ch, "y"] <- length(rows) - r + 1
    }
    coords
})()


.plot_network_v2 <- function(M, thr, file,
                             labels = TARGET_CHANNELS_V2,
                             coords = LAYOUT_V2) {
    A <- (abs(M) >= thr) * 1L
    diag(A) <- 0
    g <- igraph::graph_from_adjacency_matrix(A, mode = "directed")
    igraph::V(g)$label       <- toupper(labels)
    igraph::V(g)$size        <- 14
    igraph::V(g)$color       <- "lightblue"
    igraph::V(g)$frame.color <- "gray40"
    igraph::V(g)$label.cex   <- 0.8
    igraph::V(g)$label.color <- "darkblue"
    igraph::E(g)$arrow.size <- 0.25
    igraph::E(g)$color       <- adjustcolor("gray30", alpha.f = 0.6)

    pdf(file, width = 5, height = 4.2)
    op <- par(mar = c(0.2, 0.2, 0.2, 0.2))
    plot(g, layout = coords[labels, ],
         rescale = TRUE, margin = 0.10, asp = 0.9)
    par(op); dev.off()
}


.threshold_v2 <- function(mats, keep_frac = 0.20) {
    vals <- unlist(lapply(mats, function(M) abs(M)))
    stats::quantile(vals, 1 - keep_frac, na.rm = TRUE)
}


write_network_pdfs_v2 <- function(per_subj, ica, keep_frac = 0.20) {
    for (k in 1:3) {
        clusters <- ica$mild[[k]]$clusters
        K <- length(clusters)
        L_avg <- lapply(clusters, function(g) .cluster_average(per_subj[[k]]$left,  g))
        R_avg <- lapply(clusters, function(g) .cluster_average(per_subj[[k]]$right, g))
        thr   <- .threshold_v2(c(L_avg, R_avg), keep_frac = keep_frac)
        for (c in 1:K) {
            .plot_network_v2(L_avg[[c]], thr,
                file.path(IMG_V2, sprintf("cp%d_cluster%d_left.pdf",  k, c)))
            .plot_network_v2(R_avg[[c]], thr,
                file.path(IMG_V2, sprintf("cp%d_cluster%d_right.pdf", k, c)))
        }
        cat(sprintf("[v2 networks] CP%d -> %d clusters (thr=%.3f)\n", k, K, thr))
    }
}


write_hist_v2 <- function(fits) {
    file <- file.path(IMG_V2, "hist_cps.pdf")
    pdf(file, width = 7, height = 6)
    op <- par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
    cols <- c("lightblue", "mistyrose", "wheat")
    for (k in 1:3) {
        hist(fits[[k]]$per_sub_cp, breaks = 25, xlim = c(0, fits[[k]]$n),
             col = cols[k],
             main = sprintf("Histogram of detected change point %d", k),
             xlab = "Location", ylab = "Frequency")
        abline(v = fits[[k]]$cohort_cp, col = "red", lwd = 2, lty = 2)
    }
    par(op); dev.off()
    cat(sprintf("[v2 hist] -> %s\n", file))
}


write_hamming_pdfs_v2 <- function(H, setting = "mild") {
    palette <- colorRampPalette(c("white", "red"))(50)
    for (k in 1:3) {
        for (side in c("left", "right")) {
            M    <- H[[setting]][[k]][[paste0("H_", side)]]
            file <- file.path(IMG_V2, sprintf("cp%d_%s_hamming.pdf", k, side))
            pdf(file, width = 5, height = 4)
            pheatmap::pheatmap(M, color = palette,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               display_numbers = TRUE, number_color = "black",
                               fontsize_number = 11, legend = FALSE)
            dev.off()
        }
        cat(sprintf("[v2 hamming] CP%d  left+right\n", k))
    }
}


write_outlier_v2 <- function(intervals, fits, outlier = 13L, normal = 6L) {
    subj <- intervals[[1]]; cps <- fits[[1]]$per_sub_cp
    for (info in list(list(s = outlier, name = "subject_13.pdf"),
                      list(s = normal,  name = "subject_6.pdf"))) {
        X <- subj[[info$s]]; cp <- cps[info$s]
        pdf(file.path(IMG_V2, info$name), width = 5, height = 4)
        op <- par(mar = c(4, 4, 2, 1))
        matplot(seq_len(nrow(X)), X, type = "l", lty = 1,
                col = grDevices::rainbow(ncol(X), alpha = 0.7),
                xlab = "time", ylab = "", main = "")
        abline(v = cp, col = "red", lwd = 2)
        par(op); dev.off()
    }
    cat(sprintf("[v2 outlier] subj %d cp=%d ; subj %d cp=%d\n",
                outlier, cps[outlier], normal, cps[normal]))
}


## ============= driver =====================================================
if (sys.nframe() == 0L) {
    cat("\n>>> Paper montage v2: 19 standard 10-20 + FPZ + OZ (no mastoids)\n\n")
    ## The processed EEG (intervals.rds) ships in application/cache/, so the
    ## full analysis reproduces without the raw BDF files. Only fall back to
    ## reading raw_data/ + preprocessing if that cache is absent.
    int_cache <- file.path(CACHE_V2, "intervals.rds")
    if (file.exists(int_cache)) {
        cat(sprintf("[intervals] using shipped cache %s (raw_data/ not needed)\n",
                    int_cache))
        intervals <- readRDS(int_cache)
    } else {
        intervals <- intervals_v2(preproc_v2())
    }
    fits      <- cohort_v2(intervals)

    p <- fits[[1]]$p; n <- fits[[1]]$n
    cat("\n--- Cohort CP summary ---\n")
    for (k in 1:3) cat(sprintf("  CP%d  cohort = %d  mean = %.3f  sd = %.3f\n",
                               k, fits[[k]]$cohort_cp,
                               mean(fits[[k]]$per_sub_cp) / n,
                               sd  (fits[[k]]$per_sub_cp) / n))

    cat("\n--- ICA (gammas 200 / 375 same as main run; may need retune) ---\n")
    mild <- lapply(fits, function(f) cluster_from_sse(f$sse_list, p, n, gamma = 200))
    cons <- lapply(fits, function(f) cluster_from_sse(f$sse_list, p, n, gamma = 375))
    for (k in 1:3) cat(sprintf("  CP%d  mild #clust = %d   cons #clust = %d\n",
                               k, length(mild[[k]]$clusters), length(cons[[k]]$clusters)))
    ica <- list(mild = mild, conservative = cons)
    saveRDS(ica, file.path(CACHE_V2, "ica.rds"))

    ## Paper Tables 8-11 from the montage fits/ica, using the same writers as
    ## scripts 03/04 so the CSV format matches (these emit the v2 numbers, not
    ## the old mastoid-montage ones the standalone scripts produce).
    write_table8 (fits)
    write_table9 (ica)
    write_table10(ica, fits)
    write_table11(ica)

    per_subj <- networks_v2(intervals, fits)

    cat("\n--- Hamming under v2 ---\n")
    H <- compute_hamming_tables(per_subj, ica)
    saveRDS(H, file.path(CACHE_V2, "hamming.rds"))
    for (setting in c("mild", "conservative")) for (k in 1:3) {
        M <- round(H[[setting]][[k]]$H_avg, 3)
        cat(sprintf("\n%s  CP%d  (thr=%.3f):\n", setting, k, H[[setting]][[k]]$thr))
        print(M)
        write.csv(M, file.path(OUT_V2, sprintf("hamming_%s_cp%d.csv", setting, k)))
    }

    cat("\n--- within vs max-between Hamming per cluster ---\n")
    for (setting in c("mild", "conservative")) for (k in 1:3) {
        M <- H[[setting]][[k]]$H_avg
        d <- diag(M); offd <- M; diag(offd) <- NA
        max_o <- apply(offd, 1, max, na.rm = TRUE)
        viol  <- sum(d > max_o)
        cat(sprintf("  %-12s CP%d:  diag = [%s]  max off = [%s]  violators %d/%d\n",
                    setting, k,
                    paste(round(d,     3), collapse = ", "),
                    paste(round(max_o, 3), collapse = ", "),
                    viol, length(d)))
    }

    cat("\n--- writing publication PDFs to images/ ---\n")
    write_hist_v2        (fits)
    write_network_pdfs_v2(per_subj, ica)
    write_hamming_pdfs_v2(H, setting = "mild")
    write_outlier_v2     (intervals, fits)
}

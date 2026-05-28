### Application Section 5.2: per-cluster averaged networks (Figs 2-4).
###
### For each subject, fits sparse VAR(1) transition matrices on both sides of
### that subject's estimated change point via fista.sparse. These per-subject
### A_left, A_right matrices are cached in application/cache/per_subject_A.rds
### and consumed by both this script (network plots) and 06_hamming.R.
###
### Network plotting averages the real-valued matrices within each cluster,
### thresholds at a fixed quantile of |entries|, and renders an igraph layout.

suppressPackageStartupMessages({
    library(igraph)
})

source("auxiliary_functions.R")
source("application/04_ica_clustering.R")

FISTA_NITER <- 80
A_CACHE     <- "application/cache/per_subject_A.rds"


#' Fit sparse VAR(1) on a segment by FISTA. Returns the p x p phi.hat.
.fit_segment <- function(seg, lambda) {
    n <- nrow(seg); p <- ncol(seg)
    X <- seg[1:(n - 1), , drop = FALSE]
    Y <- seg[2:n,       , drop = FALSE]
    fit <- fista.sparse(X, Y, lambda = lambda, d = p, niter = FISTA_NITER)
    fit$phi.hat
}


#' For each subject in each interval, fit left- and right-side VAR(1) at the
#' subject's own estimated CP. Returns a length-3 list, each a list with
#' components $left and $right (each a length-M list of p x p matrices).
fit_per_subject_networks <- function(intervals, fits,
                                     lambdas = LAMBDAS, force = FALSE) {
    if (!force && file.exists(A_CACHE)) {
        cat("[networks] cached -> ", A_CACHE, "\n")
        return(readRDS(A_CACHE))
    }
    out <- vector("list", length(intervals))
    for (k in seq_along(intervals)) {
        subj <- intervals[[k]]
        cps  <- fits[[k]]$per_sub_cp
        M    <- length(subj)
        left  <- vector("list", M)
        right <- vector("list", M)
        for (i in 1:M) {
            X  <- subj[[i]]; cp <- cps[i]
            n  <- nrow(X)
            n_left  <- cp;  n_right <- n - cp
            lam_L <- if (n_left  <= n_right) lambdas[1] else lambdas[2]
            lam_R <- if (n_right <  n_left)  lambdas[1] else lambdas[2]
            left [[i]] <- .fit_segment(X[1:cp, , drop = FALSE],         lam_L)
            right[[i]] <- .fit_segment(X[(cp + 1):n, , drop = FALSE],   lam_R)
        }
        out[[k]] <- list(left = left, right = right)
        cat(sprintf("  CP%d: fit %d subjects (left + right)\n", k, M))
    }
    saveRDS(out, A_CACHE)
    invisible(out)
}


#' Average transition matrices within a cluster.
.cluster_average <- function(mat_list, members) {
    Reduce("+", mat_list[members]) / length(members)
}


#' Pick a threshold = (1 - keep_frac) quantile of |entries| over all matrices.
.threshold <- function(mats, keep_frac = 0.20) {
    vals <- unlist(lapply(mats, function(M) abs(M)))
    stats::quantile(vals, 1 - keep_frac, na.rm = TRUE)
}


#' Plot a single thresholded transition matrix as an igraph network.
#' Self-loops removed for legibility; edge alpha encodes weight.
.plot_one_network <- function(M, thr, title = "") {
    p <- nrow(M)
    g <- igraph::graph_from_adjacency_matrix(
        abs(M) >= thr, mode = "directed", diag = FALSE)
    igraph::V(g)$label   <- NA
    igraph::V(g)$size    <- 8
    igraph::V(g)$color   <- "skyblue"
    igraph::E(g)$arrow.size <- 0.15
    igraph::E(g)$color   <- adjustcolor("gray30", alpha.f = 0.45)
    plot(g, layout = igraph::layout_in_circle(g),
         vertex.frame.color = "gray40",
         vertex.label.cex   = 0.6,
         main = title, margin = c(0, 0, 0, 0))
}


#' Figs 2-4: one PDF per interval, top row = left-side cluster averages,
#' bottom row = right-side cluster averages.
write_network_figs <- function(per_subj, ica,
                               out_dir = "application/output",
                               keep_frac = 0.20) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    for (k in 1:3) {
        clusters <- ica$mild[[k]]$clusters
        K        <- length(clusters)
        L_avg <- lapply(clusters, function(g) .cluster_average(per_subj[[k]]$left,  g))
        R_avg <- lapply(clusters, function(g) .cluster_average(per_subj[[k]]$right, g))
        thr   <- .threshold(c(L_avg, R_avg), keep_frac = keep_frac)

        out <- file.path(out_dir, sprintf("fig%d_networks_cp%d.pdf", 1 + k, k))
        pdf(out, width = 3 * K, height = 6)
        op  <- par(mfrow = c(2, K), mar = c(1, 1, 2, 1))
        for (c in 1:K) .plot_one_network(L_avg[[c]], thr,
                                          sprintf("CP%d L  cluster%d", k, c))
        for (c in 1:K) .plot_one_network(R_avg[[c]], thr,
                                          sprintf("CP%d R  cluster%d", k, c))
        par(op); dev.off()
        cat(sprintf("[fig %d] -> %s   (K = %d clusters, threshold = %.4f)\n",
                    1 + k, out, K, thr))
    }
}


if (sys.nframe() == 0L) {
    intervals <- build_intervals()
    fits      <- run_cohort_cp()
    ica       <- if (file.exists("application/cache/ica.rds"))
                     readRDS("application/cache/ica.rds")
                 else
                     run_ica(fits)
    per_subj  <- fit_per_subject_networks(intervals, fits)
    write_network_figs(per_subj, ica)
}

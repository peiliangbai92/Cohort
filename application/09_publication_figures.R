### Regenerate the per-cluster, per-side network plots, per-interval Hamming
### heatmaps, the CP-histogram panel, and the subject-13 outlier panels using
### the filenames that tex/Cohort_cp_detection_manuscript/main.tex expects.

suppressPackageStartupMessages({
    library(igraph); library(pheatmap); library(grid)
})

source("application/02_build_intervals.R")
source("application/03_cohort_cp.R")
source("application/04_ica_clustering.R")
source("application/05_networks.R")
source("application/06_hamming.R")

IMG_DIR <- "tex/Cohort_cp_detection_manuscript/images"


#' One thresholded igraph network as a standalone PDF (5x5 inch).
#' Uses generous internal padding so no nodes are clipped at the bbox.
.one_network_pdf <- function(M, thr, file) {
    p <- nrow(M)
    g <- igraph::graph_from_adjacency_matrix(
        abs(M) >= thr, mode = "directed", diag = FALSE)
    igraph::V(g)$label      <- NA
    igraph::V(g)$size        <- 6
    igraph::V(g)$color       <- "skyblue"
    igraph::V(g)$frame.color <- "gray40"
    igraph::E(g)$arrow.size <- 0.15
    igraph::E(g)$color      <- adjustcolor("gray30", alpha.f = 0.45)

    pdf(file, width = 5, height = 5)
    op <- par(mar = c(0.2, 0.2, 0.2, 0.2))
    plot(g, layout = igraph::layout_in_circle(g),
         rescale = TRUE,
         xlim    = c(-1, 1), ylim = c(-1, 1),
         asp     = 1,
         margin  = 0.18)
    par(op); dev.off()
}


write_network_pdfs <- function(per_subj, ica, keep_frac = 0.20) {
    for (k in 1:3) {
        clusters <- ica$mild[[k]]$clusters
        K        <- length(clusters)
        L_avg <- lapply(clusters, function(g) .cluster_average(per_subj[[k]]$left,  g))
        R_avg <- lapply(clusters, function(g) .cluster_average(per_subj[[k]]$right, g))
        thr   <- .threshold(c(L_avg, R_avg), keep_frac = keep_frac)
        for (c in 1:K) {
            .one_network_pdf(L_avg[[c]], thr,
                file.path(IMG_DIR, sprintf("cp%d_cluster%d_left.pdf",  k, c)))
            .one_network_pdf(R_avg[[c]], thr,
                file.path(IMG_DIR, sprintf("cp%d_cluster%d_right.pdf", k, c)))
        }
        cat(sprintf("[networks] CP%d: wrote %d cluster pairs (thr=%.3f)\n", k, K, thr))
    }
}


.one_hamming_pdf <- function(M, file) {
    pdf(file, width = 5, height = 4)
    palette <- colorRampPalette(c("white", "red"))(50)
    ph <- pheatmap::pheatmap(
        M, color = palette, cluster_rows = FALSE, cluster_cols = FALSE,
        display_numbers = TRUE, number_color = "black", fontsize_number = 11,
        silent = FALSE, legend = FALSE)
    dev.off()
}


write_hamming_pdfs <- function(H, setting = "mild") {
    for (k in 1:3) {
        .one_hamming_pdf(H[[setting]][[k]]$H_left,
            file.path(IMG_DIR, sprintf("cp%d_left_hamming.pdf",  k)))
        .one_hamming_pdf(H[[setting]][[k]]$H_right,
            file.path(IMG_DIR, sprintf("cp%d_right_hamming.pdf", k)))
        cat(sprintf("[hamming] CP%d: wrote left+right heatmaps\n", k))
    }
}


write_hist_cps <- function(fits) {
    file <- file.path(IMG_DIR, "hist_cps.pdf")
    pdf(file, width = 7, height = 6)
    op   <- par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
    cols <- c("lightblue", "mistyrose", "wheat")
    for (k in 1:3) {
        hist(fits[[k]]$per_sub_cp, breaks = 25,
             xlim = c(0, fits[[k]]$n), col = cols[k],
             main = sprintf("Histogram of detected change point %d", k),
             xlab = "Location", ylab = "Frequency")
        abline(v = fits[[k]]$cohort_cp, col = "red", lwd = 2, lty = 2)
    }
    par(op); dev.off()
    cat(sprintf("[hist_cps] -> %s\n", file))
}


#' Standalone single-subject panel (used for subject_13.pdf and subject_6.pdf).
.subject_panel_pdf <- function(X, cp, file, title = NULL) {
    p <- ncol(X); n <- nrow(X)
    pdf(file, width = 5, height = 4)
    op <- par(mar = c(4, 4, 2, 1))
    matplot(seq_len(n), X, type = "l", lty = 1,
            col = grDevices::rainbow(p, alpha = 0.7),
            xlab = "time", ylab = "",
            main = if (is.null(title)) "" else title)
    abline(v = cp, col = "red", lwd = 2)
    par(op); dev.off()
}


write_outlier_panels <- function(intervals, fits,
                                 outlier = 13L, normal = NULL) {
    subj <- intervals[[1]]; cps <- fits[[1]]$per_sub_cp
    if (is.null(normal)) {
        cand <- which(abs(cps - fits[[1]]$cohort_cp) ==
                      min(abs(cps - fits[[1]]$cohort_cp)))
        normal <- if (6L %in% cand) 6L else cand[1]
    }
    .subject_panel_pdf(subj[[outlier]], cps[outlier],
                       file.path(IMG_DIR, "subject_13.pdf"))
    .subject_panel_pdf(subj[[normal]],  cps[normal],
                       file.path(IMG_DIR, "subject_6.pdf"))
    cat(sprintf("[outlier] subject_13 cp=%d ; normal subj=%d cp=%d\n",
                cps[outlier], normal, cps[normal]))
    invisible(list(outlier_cp = cps[outlier],
                   normal_subj = normal, normal_cp = cps[normal]))
}


if (sys.nframe() == 0L) {
    stopifnot(dir.exists(IMG_DIR))
    intervals <- build_intervals()
    fits      <- run_cohort_cp()
    ica       <- readRDS("application/cache/ica.rds")
    per_subj  <- fit_per_subject_networks(intervals, fits)
    H         <- readRDS("application/cache/hamming.rds")

    write_network_pdfs (per_subj, ica)
    write_hamming_pdfs (H, setting = "mild")
    write_hist_cps     (fits)
    info <- write_outlier_panels(intervals, fits)
    saveRDS(info, "application/cache/outlier_info.rds")
}

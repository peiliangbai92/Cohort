### Application Section 5.2 / 5.3: within/between Hamming distances.
###
### For each interval and each side (left/right of the per-subject CP), each
### subject's transition matrix is binarised at a global quantile-based
### threshold. Hamming distance between two p x p binary matrices is
###     sum(A != B) / (p * p).
###
### Tables 12-14 (mild) and 15-17 (conservative) report the average pairwise
### Hamming distance among subjects within and between clusters, averaged
### over the left- and right-side matrices. Figure 5 shows the per-side
### heatmaps (top row left, bottom row right).

suppressPackageStartupMessages({
    library(pheatmap)
    library(grid); library(gridExtra)
})

source("application/05_networks.R")


#' Pairwise Hamming distance between all binary matrices in a list.
.pairwise_hamming <- function(B_list) {
    M <- length(B_list); p2 <- length(B_list[[1]])
    D <- matrix(NA_real_, M, M)
    for (i in 1:M) for (j in i:M) {
        d <- sum(B_list[[i]] != B_list[[j]]) / p2
        D[i, j] <- d; D[j, i] <- d
    }
    D
}


#' Build a K x K matrix of within (diag) and between (off-diag) average
#' Hamming distances among subjects given cluster memberships.
#' Diagonal averages exclude (i, i) self pairs.
.cluster_hamming <- function(D, clusters) {
    K <- length(clusters)
    H <- matrix(NA_real_, K, K)
    for (a in 1:K) for (b in a:K) {
        ia <- clusters[[a]]; ib <- clusters[[b]]
        if (a == b) {
            if (length(ia) <= 1) { H[a, a] <- 0; next }
            vals <- D[ia, ia][upper.tri(D[ia, ia])]
        } else {
            vals <- as.numeric(D[ia, ib])
        }
        H[a, b] <- mean(vals); H[b, a] <- H[a, b]
    }
    rownames(H) <- colnames(H) <- paste0("Cluster ", 1:K)
    H
}


#' Binarise a single matrix at threshold thr on |entries|.
.binarise <- function(M, thr) (abs(M) >= thr) * 1L


compute_hamming_tables <- function(per_subj, ica, keep_frac = 0.20) {
    out <- list(mild = vector("list", 3), conservative = vector("list", 3))
    for (k in 1:3) {
        L_mats <- per_subj[[k]]$left
        R_mats <- per_subj[[k]]$right
        thr    <- .threshold(c(L_mats, R_mats), keep_frac = keep_frac)
        L_bin  <- lapply(L_mats, .binarise, thr = thr)
        R_bin  <- lapply(R_mats, .binarise, thr = thr)
        D_L    <- .pairwise_hamming(L_bin)
        D_R    <- .pairwise_hamming(R_bin)
        D_avg  <- (D_L + D_R) / 2
        for (setting in c("mild", "conservative")) {
            clusters <- if (setting == "mild")
                ica$mild[[k]]$clusters else ica$conservative[[k]]$clusters
            out[[setting]][[k]] <- list(
                H_left  = .cluster_hamming(D_L,   clusters),
                H_right = .cluster_hamming(D_R,   clusters),
                H_avg   = .cluster_hamming(D_avg, clusters),
                thr     = thr
            )
        }
    }
    out
}


write_hamming_tables <- function(H, out_dir = "application/output") {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    table_ids <- list(mild = 12:14, conservative = 15:17)
    for (setting in names(H)) {
        for (k in 1:3) {
            mat <- round(H[[setting]][[k]]$H_avg, 3)
            file <- file.path(out_dir,
                sprintf("table%d_%s_cp%d.csv", table_ids[[setting]][k], setting, k))
            write.csv(mat, file)
            cat(sprintf("Table %d (%s, CP%d, threshold=%.4f):\n",
                        table_ids[[setting]][k], setting, k,
                        H[[setting]][[k]]$thr))
            print(mat)
            cat(sprintf("[table%d] -> %s\n\n",
                        table_ids[[setting]][k], file))
        }
    }
}


#' Figure 5: 2 x 3 grid of heatmaps (rows = left/right, cols = CP1/2/3).
#' Uses the mild clustering by default; pass setting = "conservative" to swap.
write_fig5 <- function(H, setting = "mild",
                       out = "application/output/fig5_hamming_heatmaps.pdf") {
    palette <- colorRampPalette(c("white", "red"))(50)
    grobs <- list()
    for (side in c("H_left", "H_right")) {
        for (k in 1:3) {
            M <- H[[setting]][[k]][[side]]
            ph <- pheatmap::pheatmap(
                M, color = palette, cluster_rows = FALSE, cluster_cols = FALSE,
                display_numbers = TRUE, number_color = "black",
                fontsize_number = 8,
                main = sprintf("CP%d %s", k, sub("H_", "", side)),
                silent = TRUE, legend = FALSE)
            grobs[[length(grobs) + 1]] <- ph$gtable
        }
    }
    pdf(out, width = 12, height = 7)
    gridExtra::grid.arrange(grobs = grobs, nrow = 2, ncol = 3)
    dev.off()
    cat(sprintf("[fig5] -> %s\n", out))
}


if (sys.nframe() == 0L) {
    intervals <- build_intervals()
    fits      <- run_cohort_cp()
    ica       <- if (file.exists("application/cache/ica.rds"))
                     readRDS("application/cache/ica.rds")
                 else
                     run_ica(fits)
    per_subj  <- fit_per_subject_networks(intervals, fits)
    H         <- compute_hamming_tables(per_subj, ica)
    write_hamming_tables(H)
    write_fig5(H, setting = "mild")
    saveRDS(H, "application/cache/hamming.rds")
    cat("[hamming] saved -> application/cache/hamming.rds\n")
}

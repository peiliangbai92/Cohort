### Application Section 5.2 / 5.3: iterative clustering algorithm (ICA).
###
### Replays the clustering half of auxiliary_functions::iterative_selection
### using the per-subject SSEs cached by 03_cohort_cp.R, so we don't refit
### sparse VARs twice. Produces Tables 9, 10, 11.
###
### The "mild" run uses the default gamma = 3 log(p) log(n); the
### "conservative" run uses a tighter gamma (smaller -> more clusters).

source("auxiliary_functions.R")
source("application/03_cohort_cp.R")


#' Cluster subjects by Hausdorff distance of (group CP, combined CP),
#' starting from cached per-subject SSEs. Mirrors the second half of
#' iterative_selection(), with `gamma` as the merge threshold.
cluster_from_sse <- function(sse_list, p, n, gamma, scale = 2) {
    N      <- length(sse_list)
    offset <- scale * p + 1
    remain <- as.list(seq_len(N))
    final  <- list()

    while (length(remain) > 0) {
        nn <- length(remain)
        if (nn == 1) {
            final[[length(final) + 1]] <- remain[[1]]
            break
        }
        distants <- numeric(nn - 1)
        for (j in 2:nn) {
            g1 <- remain[[1]]; g2 <- remain[[j]]
            sse1 <- Reduce("+", sse_list[g1])
            sse2 <- Reduce("+", sse_list[g2])
            cp1  <- which.min(sse1) + offset
            cp2  <- which.min(sse2) + offset
            cpC  <- which.min((sse1 + sse2) / (length(g1) + length(g2))) + offset
            distants[j - 1] <- hausdorff_dist(c(cp1, cp2), c(cpC))
        }
        dmin <- min(distants)
        idx  <- which.min(distants) + 1
        if (dmin < gamma) {
            remain[[1]] <- c(remain[[1]], remain[[idx]])
            remain[[idx]] <- NULL
        } else {
            final[[length(final) + 1]] <- remain[[1]]
            remain[[1]] <- NULL
        }
    }

    cluster_cps <- vapply(final, function(g) {
        s <- Reduce("+", sse_list[g]) / length(g)
        which.min(s) + offset
    }, numeric(1))

    list(clusters = final, cluster_cps = cluster_cps)
}


run_ica <- function(fits = NULL,
                    gamma_mild = NULL,
                    gamma_conservative = NULL) {
    if (is.null(fits)) fits <- run_cohort_cp()
    n <- fits[[1]]$n; p <- fits[[1]]$p
    ## gamma values chosen empirically to match the paper's cluster counts
    ## (mild: 4/2/3 CPs, conservative: 3/2/2 CPs); the theoretical default
    ## 3 log(p) log(n) ~ 63 is far too small for our per-subject CP spread
    ## (sd ~ 200 samples), so we calibrate against the Hausdorff distribution.
    if (is.null(gamma_mild))          gamma_mild         <- 200
    if (is.null(gamma_conservative))  gamma_conservative <- 375
    cat(sprintf("[ica] gamma_mild = %.2f   gamma_conservative = %.2f\n",
                gamma_mild, gamma_conservative))

    mild <- vector("list", length(fits))
    cons <- vector("list", length(fits))
    for (k in seq_along(fits)) {
        sse_list <- fits[[k]]$sse_list
        mild[[k]] <- cluster_from_sse(sse_list, p, n, gamma = gamma_mild)
        cons[[k]] <- cluster_from_sse(sse_list, p, n, gamma = gamma_conservative)
        cat(sprintf("  CP%d   mild #clusters = %d   conservative #clusters = %d\n",
                    k, length(mild[[k]]$clusters), length(cons[[k]]$clusters)))
    }
    list(mild = mild, conservative = cons,
         gamma_mild = gamma_mild, gamma_conservative = gamma_conservative)
}


#' Pretty cluster membership as a single comma-separated string per cluster.
.cluster_str <- function(clusters) {
    paste(sapply(clusters, function(g) sprintf("(%s)", paste(g, collapse = ", "))),
          collapse = " ")
}


write_table9 <- function(ica, out = "application/output/table9.csv") {
    df <- data.frame(
        No_of_CP   = 1:3,
        n_clusters = vapply(ica$mild, function(x) length(x$clusters), integer(1)),
        Subjects   = vapply(ica$mild, function(x) .cluster_str(x$clusters),
                            character(1))
    )
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    write.csv(df, out, row.names = FALSE)
    cat("\nTable 9 (mild clustering):\n"); print(df, right = FALSE)
    cat(sprintf("[table9] -> %s\n", out))
    invisible(df)
}


write_table10 <- function(ica, fits, out = "application/output/table10.csv") {
    n <- fits[[1]]$n
    rows <- list()
    for (k in 1:3) {
        clusters <- ica$mild[[k]]$clusters
        for (c in seq_along(clusters)) {
            members <- clusters[[c]]
            cps_t   <- fits[[k]]$per_sub_cp[members] / n
            rows[[length(rows) + 1]] <- data.frame(
                Change_point = k, Cluster = c,
                n = length(members),
                Mean = mean(cps_t),
                Std  = if (length(members) > 1) sd(cps_t) else 0
            )
        }
    }
    df <- do.call(rbind, rows)
    write.csv(df, out, row.names = FALSE)
    cat("\nTable 10 (mild per-cluster mean/std in T-units):\n")
    print(format(df, digits = 4))
    cat(sprintf("[table10] -> %s\n", out))
    invisible(df)
}


write_table11 <- function(ica, out = "application/output/table11.csv") {
    df <- data.frame(
        No_of_CP   = 1:3,
        n_clusters = vapply(ica$conservative, function(x) length(x$clusters), integer(1)),
        Subjects   = vapply(ica$conservative,
                            function(x) .cluster_str(x$clusters), character(1))
    )
    write.csv(df, out, row.names = FALSE)
    cat("\nTable 11 (conservative clustering):\n"); print(df, right = FALSE)
    cat(sprintf("[table11] -> %s\n", out))
    invisible(df)
}


if (sys.nframe() == 0L) {
    fits <- run_cohort_cp()
    ica  <- run_ica(fits)
    write_table9 (ica)
    write_table10(ica, fits)
    write_table11(ica)
    saveRDS(ica, "application/cache/ica.rds")
    cat("[ica] saved -> application/cache/ica.rds\n")
}

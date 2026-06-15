### Aggregate Group (G1-G4) clustering scenario results into the format of
### Table 5 in the manuscript:
###   scenario, truth (g, tau_j*/T), \hat g (mean #clusters), per-CP mean/sd,
###   global accuracy.

repo_dir <- getwd()   # run from the repo root
setwd(repo_dir)

scenario_meta <- list(
    G1 = list(T = 200, M = 20, group_sizes = c(5, 15),       taus = c(67, 133)),
    G2 = list(T = 200, M = 20, group_sizes = c(10, 10),      taus = c(67, 133)),
    G3 = list(T = 400, M = 25, group_sizes = c(5, 8, 12),    taus = c(100, 200, 300)),
    G4 = list(T = 400, M = 20, group_sizes = c(3, 7, 10),    taus = c(50, 200, 320))
)

global_accuracy <- function(cluster_assignments, group_sizes) {
    ## cluster_assignments: integer vector of length M with cluster labels
    ## group_sizes: vector of true group sizes; subjects 1..g1 belong to true
    ## group 1, then g1+1..g1+g2 to group 2, etc.
    M <- sum(group_sizes)
    true_label <- rep(seq_along(group_sizes), group_sizes)
    est_label <- cluster_assignments
    if (length(est_label) != M) return(NA)

    ## Build all clusters and find best 1-1 assignment from estimated to true
    ## (Hungarian-like by frequency; for small cardinality we enumerate).
    n_est <- length(unique(est_label))
    n_true <- length(group_sizes)
    if (n_est == 0) return(0)

    ## Confusion matrix
    conf <- table(factor(est_label, levels = sort(unique(est_label))),
                  factor(true_label, levels = seq_along(group_sizes)))

    ## Greedy assignment by maximum row entry (works when n_est >= n_true)
    perm <- integer(nrow(conf))
    used_cols <- integer(0)
    for (r in seq_len(nrow(conf))) {
        avail_cols <- setdiff(seq_len(ncol(conf)), used_cols)
        if (length(avail_cols) == 0) {
            perm[r] <- NA
        } else {
            c_best <- avail_cols[which.max(conf[r, avail_cols])]
            perm[r] <- c_best
            used_cols <- c(used_cols, c_best)
        }
    }

    ## Compute TP, TN, FP, FN over all (subject, cluster) pairs
    g <- max(n_true, n_est)
    ttp <- ttn <- tfp <- tfn <- 0
    for (i in seq_len(M)) {
        ti <- true_label[i]
        ei <- est_label[i]
        for (k in seq_len(g)) {
            in_true <- (ti == k)
            mapped_k <- if (k <= length(perm)) which(perm == k) else integer(0)
            in_est <- (length(mapped_k) > 0 && ei == mapped_k[1])
            if (in_true && in_est)        ttp <- ttp + 1
            else if (!in_true && !in_est) ttn <- ttn + 1
            else if (!in_true && in_est)  tfp <- tfp + 1
            else                          tfn <- tfn + 1
        }
    }
    (ttp + ttn) / (ttp + ttn + tfp + tfn)
}

summarise_G <- function(name, meta) {
    fcp_path  <- sprintf("simulations/%s/final_est_cp.RData", name)
    clu_path  <- sprintf("simulations/%s/final_cluster.RData", name)
    if (!file.exists(fcp_path) || !file.exists(clu_path)) {
        return(data.frame(scenario = name))
    }
    env <- new.env()
    load(fcp_path, envir = env)
    load(clu_path, envir = env)
    cp_list   <- env$cp_result        # list(niter) of estimated CPs vectors
    clu_list  <- env$cluster_result   # list(niter) of lists of subject indices

    T <- meta$T
    taus <- meta$taus
    M <- meta$M
    n_true <- length(taus)

    ## Convert clusters (list of subject-index vectors) to assignment vector
    cluster_to_assign <- function(clu) {
        if (is.null(clu) || length(clu) == 0) return(NULL)
        ass <- integer(M)
        for (k in seq_along(clu)) {
            members <- clu[[k]]
            if (!is.null(members) && length(members) > 0)
                ass[members] <- k
        }
        if (any(ass == 0)) return(NULL)
        ass
    }

    g_hats <- integer(0)
    accs   <- numeric(0)
    cp_per_cluster <- vector("list", n_true)   # one bucket per true CP
    for (j in seq_len(n_true)) cp_per_cluster[[j]] <- numeric(0)

    for (e in seq_along(cp_list)) {
        cps_e <- cp_list[[e]]
        ass_e <- cluster_to_assign(clu_list[[e]])
        if (is.null(cps_e)) next

        g_hats <- c(g_hats, length(cps_e))

        ## Bucket each estimated CP into the closest true CP
        cps_sorted <- sort(cps_e)
        for (cp_val in cps_sorted) {
            j <- which.min(abs(cp_val - taus))
            cp_per_cluster[[j]] <- c(cp_per_cluster[[j]], cp_val)
        }

        if (!is.null(ass_e)) {
            accs <- c(accs, global_accuracy(ass_e, meta$group_sizes))
        }
    }

    truth_str <- sprintf("(%d, %s)", n_true,
                         paste(sprintf("%.3f", taus / T), collapse = ", "))
    mean_str  <- paste(sprintf("%.4f",
                               sapply(cp_per_cluster, function(x)
                                   if (length(x)) mean(x) / T else NA)),
                       collapse = ", ")
    sd_str    <- paste(sprintf("%.4f",
                               sapply(cp_per_cluster, function(x)
                                   if (length(x) >= 2) sd(x) / T else NA)),
                       collapse = ", ")
    data.frame(scenario = name, truth = truth_str, g_hat = mean(g_hats),
               mean = mean_str, sd = sd_str,
               accuracy = round(mean(accs), 3))
}

g_summary <- do.call(rbind, lapply(names(scenario_meta),
                                    function(n) summarise_G(n, scenario_meta[[n]])))
print(g_summary, row.names = FALSE)
saveRDS(g_summary, file = "logs/summary_G.rds")
write.csv(g_summary, file = "logs/summary_G.csv", row.names = FALSE)
cat("\nSaved to logs/summary_G.{rds,csv}\n")

### Aggregate comparison results (SBS/BSS/fvarseg) into Table 7 / 8 format.
### Reads whatever RData is available in comparison/.

repo_dir <- "/Users/peiliangbai/Documents/GitHub/Cohort"
setwd(repo_dir)

## ---- helpers ------------------------------------------------------------

## Global accuracy matching Section "Performance evaluation" of the manuscript.
## true_label[i] is the true group of subject i.  est_label[i] is the cluster
## label assigned by the method (0 if the subject had no CP detected).
global_accuracy <- function(est_label, true_label) {
    M <- length(true_label)
    if (length(est_label) != M) return(NA)
    if (all(est_label == 0))     return(0)
    n_true <- length(unique(true_label))

    ## confusion table; rows = estimated cluster labels (>0), columns = truth
    nz <- est_label != 0
    if (!any(nz)) return(0)
    conf <- table(factor(est_label[nz], levels = sort(unique(est_label[nz]))),
                  factor(true_label[nz], levels = seq_len(n_true)))

    ## greedy 1-to-1 mapping from estimated cluster to true group
    perm <- integer(nrow(conf))
    used <- integer(0)
    for (r in seq_len(nrow(conf))) {
        avail <- setdiff(seq_len(ncol(conf)), used)
        if (length(avail) == 0) { perm[r] <- NA; next }
        c_best <- avail[which.max(conf[r, avail])]
        perm[r] <- c_best
        used <- c(used, c_best)
    }

    g <- max(n_true, nrow(conf))
    ttp <- ttn <- tfp <- tfn <- 0
    for (i in seq_len(M)) {
        ti <- true_label[i]
        ei <- est_label[i]
        for (k in seq_len(g)) {
            in_true <- (ti == k)
            mapped_k <- which(perm == k)
            in_est  <- (length(mapped_k) > 0 && ei == as.integer(rownames(conf))[mapped_k[1]])
            if (in_true && in_est)        ttp <- ttp + 1
            else if (!in_true && !in_est) ttn <- ttn + 1
            else if (!in_true && in_est)  tfp <- tfp + 1
            else                          tfn <- tfn + 1
        }
    }
    (ttp + ttn) / (ttp + ttn + tfp + tfn)
}


summarise_one <- function(rdata_path, method_name) {
    if (!file.exists(rdata_path)) {
        return(data.frame(method = method_name, niter_used = 0,
                          truth = NA_character_, g_hat = NA_real_,
                          mean = NA_character_, sd = NA_character_,
                          accuracy = NA_real_))
    }
    env <- new.env()
    load(rdata_path, envir = env)
    all_cps     <- env$all_cps
    all_clust   <- env$all_clust
    all_subj_id <- env$all_subj_id
    tau_true    <- env$tau_true
    group_sizes <- env$group_sizes
    T <- env$T; M <- env$M

    true_label <- rep(seq_along(group_sizes), group_sizes)
    n_true <- length(tau_true)
    niter <- length(all_cps)
    ## consider only epochs that ran (some entries might be NULL if partial)
    valid <- which(vapply(all_cps, function(x) !is.null(x), logical(1)))
    if (length(valid) == 0) {
        return(data.frame(method = method_name, niter_used = 0,
                          truth = NA_character_, g_hat = NA_real_,
                          mean = NA_character_, sd = NA_character_,
                          accuracy = NA_real_))
    }

    cp_per_truth <- vector("list", n_true)
    for (j in seq_len(n_true)) cp_per_truth[[j]] <- numeric(0)

    g_hats <- integer(0)
    accs   <- numeric(0)

    for (e in valid) {
        cps_e   <- all_cps[[e]]
        clust_e <- all_clust[[e]]
        sub_id  <- all_subj_id[[e]]

        ## number of nonempty clusters
        if (length(clust_e) == 0) {
            g_hats <- c(g_hats, 0)
        } else {
            g_hats <- c(g_hats, length(setdiff(unique(clust_e), 0)))
        }

        ## bucket each cp by nearest truth
        if (length(cps_e) > 0) {
            for (cp_val in cps_e) {
                j <- which.min(abs(cp_val - tau_true))
                cp_per_truth[[j]] <- c(cp_per_truth[[j]], cp_val)
            }
        }

        ## per-subject assignment for accuracy: take majority cluster label
        ## of each subject's reported cps
        ass <- integer(M)
        if (length(cps_e) > 0) {
            for (i in seq_len(M)) {
                idx <- which(sub_id == i)
                if (length(idx) > 0) {
                    ass[i] <- as.integer(names(sort(table(clust_e[idx]),
                                                    decreasing = TRUE))[1])
                }
            }
        }
        accs <- c(accs, global_accuracy(ass, true_label))
    }

    truth_str <- sprintf("(%d, %s)", n_true,
                         paste(sprintf("%.3f", tau_true / T), collapse = ", "))
    mean_str  <- paste(sprintf("%.4f",
                               sapply(cp_per_truth, function(x)
                                   if (length(x)) mean(x) / T else NA)),
                       collapse = ", ")
    sd_str    <- paste(sprintf("%.4f",
                               sapply(cp_per_truth, function(x)
                                   if (length(x) >= 2) sd(x) / T else NA)),
                       collapse = ", ")
    data.frame(method     = method_name,
               niter_used = length(valid),
               truth      = truth_str,
               g_hat      = round(mean(g_hats), 2),
               mean       = mean_str,
               sd         = sd_str,
               accuracy   = round(mean(accs, na.rm = TRUE), 3))
}


## ---- run ----------------------------------------------------------------
res <- rbind(
    summarise_one("comparison/compare_SBS.RData",     "SBS+kmeans"),
    summarise_one("comparison/compare_BSS.RData",     "BSS+kmeans"),
    summarise_one("comparison/compare_fvarseg.RData", "FVARSeg+kmeans")
)
print(res, row.names = FALSE)
saveRDS(res, file = "logs/comparison_summary.rds")
write.csv(res, file = "logs/comparison_summary.csv", row.names = FALSE)
cat("\nSaved to logs/comparison_summary.{rds,csv}\n")

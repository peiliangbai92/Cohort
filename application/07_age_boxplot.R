### Application Section 5.3 (Figure 6): age boxplot + permutation t-test
### across the two largest clusters of the conservative CP1 partition.
###
### Subject ages live in application/data/subject_ages.csv with columns
### "subject" (1..22) and "age". This file is NOT part of the BDF release;
### it lives in the rsed2017 participant metadata. If it is missing the
### script writes a stub CSV and aborts with a clear message; fill it in and
### re-run.

source("application/04_ica_clustering.R")

AGES_CSV   <- "application/data/subject_ages.csv"
N_PERM     <- 10000


.load_ages <- function(path = AGES_CSV) {
    if (!file.exists(path)) {
        dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
        stub <- data.frame(subject = 1:22, age = NA_real_)
        write.csv(stub, path, row.names = FALSE)
        stop(sprintf(
            "Missing ages file. A stub has been written to %s.\nFill the 'age' column from the rsed2017 metadata and re-run.",
            path))
    }
    df <- read.csv(path)
    if (anyNA(df$age)) stop("Some ages are NA in ", path)
    df
}


#' Two-sample permutation t-test on age means (two-sided).
permutation_t_test <- function(a, b, n_perm = N_PERM, seed = 1L) {
    set.seed(seed)
    pooled <- c(a, b); na <- length(a)
    obs <- (mean(a) - mean(b)) /
        sqrt(var(a) / na + var(b) / (length(pooled) - na))
    null_stats <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
        idx  <- sample.int(length(pooled), na)
        a_p  <- pooled[idx]; b_p <- pooled[-idx]
        null_stats[i] <- (mean(a_p) - mean(b_p)) /
            sqrt(var(a_p) / na + var(b_p) / (length(pooled) - na))
    }
    list(t = obs, p.value = mean(abs(null_stats) >= abs(obs)))
}


#' Use the conservative partition of CP1, keep the two clusters with >= 5
#' members (matches the paper's Section 5.3 description), draw the boxplot,
#' run the permutation t-test, and return the test object.
write_fig6 <- function(ica = NULL, ages_df = NULL,
                       out = "application/output/fig6_age_boxplot.pdf",
                       min_cluster_size = 5) {
    if (is.null(ica))     ica     <- readRDS("application/cache/ica.rds")
    if (is.null(ages_df)) ages_df <- .load_ages()

    clusters <- ica$conservative[[1]]$clusters
    big      <- clusters[lengths(clusters) >= min_cluster_size]
    if (length(big) != 2)
        warning(sprintf("expected 2 big clusters, got %d", length(big)))

    grp_ages <- lapply(big, function(g) ages_df$age[match(g, ages_df$subject)])
    test     <- permutation_t_test(grp_ages[[1]], grp_ages[[2]])

    pdf(out, width = 5, height = 5)
    boxplot(grp_ages, names = paste0("cluster", seq_along(grp_ages)),
            ylab = "age")
    dev.off()
    cat(sprintf("[fig6] -> %s\n", out))
    cat(sprintf("  cluster sizes: %s\n",
                paste(lengths(big), collapse = ", ")))
    cat(sprintf("  cluster means: %s\n",
                paste(round(sapply(grp_ages, mean), 2), collapse = ", ")))
    cat(sprintf("  perm t = %.4f   p = %.4f\n", test$t, test$p.value))
    invisible(test)
}


if (sys.nframe() == 0L) write_fig6()

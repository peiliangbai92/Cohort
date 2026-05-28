### Application Section 5.2: cohort change-point detection on the 3 intervals.
###
### Calls auxiliary_functions::L() on each subject to get the per-subject
### sparse-VAR(1) SSE curve, then:
###   * cohort CP   = argmin over t of  sum_i SSE_i(t)       (Eq. 3 of the paper)
###   * per-subject = argmin over t of  SSE_i(t)              (Section 2.2.1)
### Produces Table 8 (mean/std of per-subject CPs + cohort estimate) and Fig 1
### (3-panel histogram of detected per-subject CPs).
###
### Caches the per-subject SSEs in application/cache/sse_lists.rds so the
### downstream ICA script can reuse them without refitting.

source("auxiliary_functions.R")
source("application/02_build_intervals.R")

LAMBDAS    <- c(1, 1)   # lambdas[1] = shorter-side L1; lambdas[2] = longer-side
SCALE      <- 2         # L() skip parameter; CP search range = [scale*p+1, T-(scale*p+1)]
SSE_CACHE  <- "application/cache/sse_lists.rds"


#' For one length-M list of subject matrices, return:
#'   sse_list   : list of length M of per-subject SSE curves (length = T - 2*(scale*p+1) + 1)
#'   per_sub_cp : vector of length M giving argmin per subject in absolute T units
#'   cohort_cp  : argmin of sum_i SSE_i, in absolute T units
fit_one_interval <- function(subjects, lambdas = LAMBDAS, scale = SCALE) {
    M  <- length(subjects)
    p  <- ncol(subjects[[1]])
    n  <- nrow(subjects[[1]])
    offset <- scale * p + 1

    sse_list <- vector("list", M)
    per_sub  <- numeric(M)
    for (i in seq_len(M)) {
        t0  <- proc.time()
        sse <- L(subjects[[i]], lambda = lambdas, scale = scale)
        dt  <- (proc.time() - t0)[3]
        sse_list[[i]] <- sse
        per_sub[i]    <- which.min(sse) + offset
        cat(sprintf("  subject %2d: cp = %4d   (L took %.1fs)\n", i, per_sub[i], dt))
    }
    combined  <- Reduce("+", sse_list)
    cohort_cp <- which.min(combined) + offset
    list(sse_list = sse_list, per_sub_cp = per_sub,
         cohort_cp = cohort_cp, combined_sse = combined,
         offset = offset, n = n, p = p)
}


run_cohort_cp <- function(intervals = NULL, force = FALSE) {
    if (!force && file.exists(SSE_CACHE)) {
        cat("[cohort_cp] cached -> ", SSE_CACHE, "\n")
        return(readRDS(SSE_CACHE))
    }
    if (is.null(intervals)) intervals <- build_intervals()

    fits <- vector("list", length(intervals))
    for (k in seq_along(intervals)) {
        cat(sprintf("\n=== Interval %d (M = %d subjects) ===\n",
                    k, length(intervals[[k]])))
        fits[[k]] <- fit_one_interval(intervals[[k]])
        cat(sprintf("  cohort cp = %d   per-subject mean = %.1f   sd = %.1f\n",
                    fits[[k]]$cohort_cp,
                    mean(fits[[k]]$per_sub_cp),
                    sd(fits[[k]]$per_sub_cp)))
    }
    saveRDS(fits, SSE_CACHE)
    invisible(fits)
}


#' Table 8: mean / sd of per-subject CPs (in T-units, i.e. divided by n),
#' plus the cohort estimate (also in T-units).
write_table8 <- function(fits, out = "application/output/table8.csv") {
    n <- fits[[1]]$n
    df <- data.frame(
        Statistic = c("Mean", "Std", "Cohort est."),
        Changepoint_1 = c(mean(fits[[1]]$per_sub_cp) / n,
                          sd(fits[[1]]$per_sub_cp)   / n,
                          fits[[1]]$cohort_cp        / n),
        Changepoint_2 = c(mean(fits[[2]]$per_sub_cp) / n,
                          sd(fits[[2]]$per_sub_cp)   / n,
                          fits[[2]]$cohort_cp        / n),
        Changepoint_3 = c(mean(fits[[3]]$per_sub_cp) / n,
                          sd(fits[[3]]$per_sub_cp)   / n,
                          fits[[3]]$cohort_cp        / n)
    )
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    write.csv(df, out, row.names = FALSE)
    cat("\nTable 8 (per-subject mean/std and cohort est. in T-units):\n")
    print(format(df, digits = 4))
    cat(sprintf("[table8] -> %s\n", out))
    invisible(df)
}


#' Figure 1: 3-panel histogram of detected per-subject change points.
write_fig1 <- function(fits, out = "application/output/fig1_histograms.pdf",
                       breaks = 25) {
    pdf(out, width = 7, height = 6)
    op <- par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
    cols <- c("lightblue", "mistyrose", "wheat")
    for (k in 1:3) {
        hist(fits[[k]]$per_sub_cp,
             breaks = breaks,
             xlim   = c(0, fits[[k]]$n),
             col    = cols[k],
             main   = sprintf("Histogram of detected change point %d", k),
             xlab   = "Location",
             ylab   = "Frequency")
        abline(v = fits[[k]]$cohort_cp, col = "red", lwd = 2, lty = 2)
    }
    par(op); dev.off()
    cat(sprintf("[fig1] -> %s\n", out))
}


if (sys.nframe() == 0L) {
    args  <- commandArgs(trailingOnly = TRUE)
    force <- "--force" %in% args
    fits  <- run_cohort_cp(force = force)
    write_table8(fits)
    write_fig1(fits)
}

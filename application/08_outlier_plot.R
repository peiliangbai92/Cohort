### Application Section 5.4 (Figure 7): time-series view of an outlier subject
### (subject 13 in the paper, detected CP around 230) next to a typical subject
### whose detected CP sits near the cohort estimate (~510).

source("application/03_cohort_cp.R")


#' Pick the typical subject for interval 1 = subject whose per-subject CP is
#' closest to the cohort estimate, excluding the outlier.
pick_typical <- function(per_sub_cp, cohort_cp, exclude = integer()) {
    diffs <- abs(per_sub_cp - cohort_cp)
    diffs[exclude] <- Inf
    which.min(diffs)
}


plot_outlier_panel <- function(X, cp, title) {
    p <- ncol(X); n <- nrow(X)
    matplot(seq_len(n), X, type = "l", lty = 1,
            col = grDevices::rainbow(p, alpha = 0.7),
            xlab = "time", ylab = "", main = title)
    abline(v = cp, col = "red", lwd = 2)
}


write_fig7 <- function(intervals       = NULL,
                       fits            = NULL,
                       outlier_subject = 13,
                       out             = "application/output/fig7_outlier.pdf") {
    if (is.null(intervals)) intervals <- build_intervals()
    if (is.null(fits))      fits      <- run_cohort_cp()
    subj    <- intervals[[1]]
    cps     <- fits[[1]]$per_sub_cp
    typical <- pick_typical(cps, fits[[1]]$cohort_cp, exclude = outlier_subject)

    pdf(out, width = 10, height = 4)
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    plot_outlier_panel(subj[[outlier_subject]], cps[outlier_subject],
                       sprintf("subj %d (outlier)  cp=%d",
                               outlier_subject, cps[outlier_subject]))
    plot_outlier_panel(subj[[typical]], cps[typical],
                       sprintf("subj %d (typical)  cp=%d",
                               typical, cps[typical]))
    par(op); dev.off()
    cat(sprintf("[fig7] -> %s   (outlier subj %d cp=%d ; typical subj %d cp=%d)\n",
                out, outlier_subject, cps[outlier_subject],
                typical, cps[typical]))
}


if (sys.nframe() == 0L) write_fig7()

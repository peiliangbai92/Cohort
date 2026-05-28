### Application Section 5.1 (cont.): build the 3 x 22 cohort intervals.
###
### For each of the three empirical change points (at samples 16320, 32640,
### 48960 in the original 256 Hz time base), cut +/- 4000 samples and
### subsample every 8 to obtain a T = 1000 window with the empirical CP at
### sample 500. Returns a length-3 list (one per CP), each a length-22 list
### of T x p matrices.

source("application/01_preprocess_bdf.R")

EMPIRICAL_CPS <- c(16320, 32640, 48960)
HALF_WIN      <- 4000          # +/- around the CP (in 256 Hz samples)
SUBSAMPLE     <- 8             # keep every 8th sample -> T = 1000


#' Window one preprocessed subject around an empirical CP.
window_subject <- function(X, cp, half = HALF_WIN, step = SUBSAMPLE) {
    lo <- cp - half + 1
    hi <- cp + half
    stopifnot(lo >= 1, hi <= nrow(X))
    seg <- X[lo:hi, , drop = FALSE]
    idx <- seq(1, nrow(seg), by = step)
    seg[idx, , drop = FALSE]
}


build_intervals <- function(preproc       = NULL,
                            cache_dir     = "application/cache",
                            intervals_file = "application/cache/intervals.rds",
                            force          = FALSE) {
    if (!force && file.exists(intervals_file)) {
        cat("[intervals] cached -> ", intervals_file, "\n")
        return(readRDS(intervals_file))
    }
    if (is.null(preproc)) preproc <- preprocess_all(cache_dir = cache_dir)
    intervals <- vector("list", length(EMPIRICAL_CPS))
    for (k in seq_along(EMPIRICAL_CPS)) {
        cp <- EMPIRICAL_CPS[k]
        subj <- lapply(preproc, window_subject, cp = cp)
        intervals[[k]] <- subj
        cat(sprintf("[intervals] CP%d @ t=%d  -> %d subjects x (%d x %d)\n",
                    k, cp, length(subj), nrow(subj[[1]]), ncol(subj[[1]])))
    }
    saveRDS(intervals, intervals_file)
    invisible(intervals)
}


if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    force <- "--force" %in% args
    build_intervals(force = force)
}

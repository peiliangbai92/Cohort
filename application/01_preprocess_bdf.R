### Application Section 5.1: BDF -> 21-channel alpha-band-filtered, detrended.
###
### Reads each of the 22 raw BDF files in raw_data/, keeps the standard 10-20
### montage (19 EEG + 2 mastoids = 21 channels; NAS is intentionally excluded
### per the paper), takes the first 4 minutes (T = 61440 samples at 256 Hz),
### applies a robust iterative polynomial detrending (de Cheveigne 2018 style),
### and band-passes to 8-13 Hz (alpha) with a zero-phase Butterworth filter.
###
### Caches the preprocessed (T x p) matrix per subject under application/cache/.
### Loaded downstream by 02_build_intervals.R.

suppressPackageStartupMessages({
    library(edfReader)
    library(signal)
})

FS               <- 256                     # sampling rate (Hz)
DURATION_S       <- 4 * 60                  # first 4 minutes
T_FULL           <- FS * DURATION_S         # 61440 samples
ALPHA_BAND       <- c(8, 13)                # Hz
DETREND_ORDER    <- 3
DETREND_MAX_ITER <- 4
DETREND_THRESH   <- 3                       # MAD multiples

## Standard 10-20 EEG montage (19 channels) plus the two mastoids = 21 channels.
## The paper drops the NAS reference channel; we never include it.
TARGET_CHANNELS <- c(
    "Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8",
    "T7",  "C3", "Cz", "C4", "T8",
    "P7",  "P3", "Pz", "P4", "P8",
    "O1",  "O2",
    "M1",  "M2"
)
stopifnot(length(TARGET_CHANNELS) == 21)


#' Robust iterative polynomial detrend (de Cheveigne 2018 style).
#'
#' Fits a polynomial of `order` to each channel by weighted least squares,
#' down-weights samples whose residual exceeds `thresh * MAD` to zero, and
#' iterates. Returns the channel minus the final fitted trend.
robust_detrend <- function(x, order = DETREND_ORDER,
                           max_iter = DETREND_MAX_ITER,
                           thresh = DETREND_THRESH) {
    n <- length(x)
    t <- seq_len(n) / n
    B <- cbind(1, stats::poly(t, degree = order, raw = FALSE))
    w <- rep(1, n)
    for (it in seq_len(max_iter)) {
        Bw <- B * w
        beta <- tryCatch(qr.solve(crossprod(Bw, B), crossprod(Bw, w * x)),
                         error = function(e) NULL)
        if (is.null(beta)) break
        fitted <- as.numeric(B %*% beta)
        res <- x - fitted
        mad_r <- stats::mad(res, constant = 1.4826)
        if (mad_r == 0) break
        w_new <- as.numeric(abs(res) <= thresh * mad_r)
        if (all(w_new == w)) break
        w <- w_new
    }
    x - fitted
}


#' Zero-phase Butterworth band-pass (8-13 Hz alpha).
alpha_bandpass <- function(x, fs = FS, band = ALPHA_BAND, order = 4) {
    nyq <- fs / 2
    bf  <- signal::butter(order, band / nyq, type = "pass")
    as.numeric(signal::filtfilt(bf, x))
}


#' Read one BDF, return T_FULL x length(TARGET_CHANNELS) matrix of
#' preprocessed (detrended + alpha-band-pass) signals.
preprocess_one_bdf <- function(path,
                               channels = TARGET_CHANNELS,
                               n_keep = T_FULL) {
    hdr <- edfReader::readEdfHeader(path)
    labels <- hdr$sHeaders$label
    miss <- setdiff(channels, labels)
    if (length(miss) > 0) {
        stop(sprintf("%s: missing channels %s",
                     basename(path), paste(miss, collapse = ", ")))
    }
    sig <- edfReader::readEdfSignals(hdr, signals = channels, simplify = FALSE)

    out <- matrix(0, nrow = n_keep, ncol = length(channels))
    colnames(out) <- channels
    for (j in seq_along(channels)) {
        v <- sig[[channels[j]]]$signal
        if (length(v) < n_keep) {
            stop(sprintf("%s/%s: only %d samples, need %d",
                         basename(path), channels[j], length(v), n_keep))
        }
        v <- v[seq_len(n_keep)]
        v <- robust_detrend(v)
        v <- alpha_bandpass(v)
        out[, j] <- v
    }
    out
}


preprocess_all <- function(raw_dir   = "raw_data",
                           cache_dir = "application/cache",
                           force     = FALSE) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    out <- vector("list", 22)
    for (sid in 1:22) {
        cache_file <- file.path(cache_dir, sprintf("preproc_S%02d.rds", sid))
        if (!force && file.exists(cache_file)) {
            out[[sid]] <- readRDS(cache_file)
            cat(sprintf("[%2d/22] cached  -> %s\n", sid, cache_file))
            next
        }
        bdf <- file.path(raw_dir, sprintf("EEG_Cat_Study4_Resting_S%d.bdf", sid))
        t0 <- proc.time()
        X  <- preprocess_one_bdf(bdf)
        dt <- (proc.time() - t0)[3]
        saveRDS(X, cache_file)
        out[[sid]] <- X
        cat(sprintf("[%2d/22] preprocess %s  (%.1fs, %d x %d)\n",
                    sid, basename(bdf), dt, nrow(X), ncol(X)))
    }
    invisible(out)
}


if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    force <- "--force" %in% args
    preprocess_all(force = force)
}

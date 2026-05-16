### Detection strategies for the cohort single change point problem.
### Implements the two estimators analyzed in Theorem 1 (max-lag) and
### Theorem 1' (min-lag) of the manuscript.
###
### Both functions assume each subject's time series is a (T x p) matrix and
### that all subjects share the cohort change point structure described in
### Section 2.1 of the manuscript. The returned estimator is the loss-function
### argmin, i.e. cp = argmin_{t in [skip, T-skip]} sum_i (left SSE + right SSE).
###
### Author: cohort change point project.

suppressPackageStartupMessages({
    library("sparsevar")
})


#' Max-lag detection strategy (Theorem 1).
#'
#' Expands all subjects to the common maximum lag `d = max_i d_i`, fits a
#' LASSO VAR(d) on each side of every candidate change point, and returns the
#' argmin of the cohort-summed SSE.
#'
#' @param subjects list of T x p matrices.
#' @param d integer, the maximum lag across subjects.
#' @param skip number of boundary time points to skip when searching.
#' @param nlambda CV grid size for `sparsevar::fitVAR`.
#' @param nfolds CV folds for `sparsevar::fitVAR`.
#' @return list with components est.cp, left_mats, right_mats, obj_vals, runtime.
single.cp.detect.maxlag <- function(subjects, d, skip = 5,
                                    nlambda = 3, nfolds = 3) {
    M <- length(subjects)
    n <- nrow(subjects[[1]])
    p <- ncol(subjects[[1]])
    stopifnot(d >= 1, skip >= d + 1, n - skip >= skip + 1)

    n_grid <- n - 2 * skip + 1
    obj_vals <- numeric(n_grid)
    start <- proc.time()
    for (k in seq_len(n_grid)) {
        t <- skip + k - 1
        sse <- 0
        for (i in 1:M) {
            X <- subjects[[i]]
            fit_left  <- fitVAR(X[1:t, , drop = FALSE], p = d,
                                nlambda = nlambda, nfolds = nfolds,
                                threshold = TRUE)
            fit_right <- fitVAR(X[(t + 1):n, , drop = FALSE], p = d,
                                nlambda = nlambda, nfolds = nfolds,
                                threshold = TRUE)
            sse <- sse + sum(fit_left$residuals^2) + sum(fit_right$residuals^2)
        }
        obj_vals[k] <- sse
    }
    elapsed <- (proc.time() - start)[3]
    cp <- skip + which.min(obj_vals) - 1

    left_mats  <- vector("list", M)
    right_mats <- vector("list", M)
    for (i in 1:M) {
        X <- subjects[[i]]
        fit_L <- fitVAR(X[1:cp, , drop = FALSE], p = d, threshold = TRUE)
        fit_R <- fitVAR(X[(cp + 1):n, , drop = FALSE], p = d, threshold = TRUE)
        left_mats[[i]]  <- fit_L$A
        right_mats[[i]] <- fit_R$A
    }

    list(est.cp = cp, left_mats = left_mats, right_mats = right_mats,
         obj_vals = obj_vals, runtime = unname(elapsed))
}


#' Truncated SSE used by `single.cp.detect.minlag`.
#'
#' Given a segment X (n0 x p), the full estimator A_list of length d_i, and a
#' truncation lag d <= d_i, returns the SSE of the truncated VAR(d) predictor
#' \hat X_t = sum_{q=1}^d A_list[[q]] X_{t-q}.
.truncated_sse <- function(X, A_list, d) {
    d_i <- length(A_list)
    stopifnot(d >= 1, d <= d_i)
    n0 <- nrow(X)
    if (n0 <= d_i) return(0)
    p <- ncol(X)
    Xhat <- matrix(0, nrow = n0 - d_i, ncol = p)
    for (q in seq_len(d)) {
        Xlag <- X[(d_i + 1 - q):(n0 - q), , drop = FALSE]
        Xhat <- Xhat + Xlag %*% t(A_list[[q]])
    }
    res <- X[(d_i + 1):n0, , drop = FALSE] - Xhat
    sum(res^2)
}


#' Min-lag detection strategy (Theorem 1').
#'
#' For each subject, fits a LASSO VAR(d_i) on both sides of every candidate
#' time point, then evaluates the cohort loss using only the truncated
#' predictor that retains the first d lags. Assumes the jump-size beyond lag
#' d vanishes, as encoded in Assumption (H2) of the manuscript.
#'
#' @param subjects list of T x p matrices.
#' @param d integer, the truncation lag (number of lags kept in the loss).
#' @param lags integer vector of length M with lags[i] = d_i >= d.
#' @param skip boundary skip.
#' @return list with components est.cp, obj_vals, runtime.
single.cp.detect.minlag <- function(subjects, d, lags, skip = 5,
                                    nlambda = 3, nfolds = 3) {
    M <- length(subjects)
    n <- nrow(subjects[[1]])
    p <- ncol(subjects[[1]])
    stopifnot(length(lags) == M, all(lags >= d), d >= 1)
    stopifnot(skip >= max(lags) + 1, n - skip >= skip + 1)

    n_grid <- n - 2 * skip + 1
    obj_vals <- numeric(n_grid)
    start <- proc.time()
    for (k in seq_len(n_grid)) {
        t <- skip + k - 1
        sse <- 0
        for (i in 1:M) {
            X <- subjects[[i]]
            d_i <- lags[i]
            fit_L <- fitVAR(X[1:t, , drop = FALSE],     p = d_i,
                            nlambda = nlambda, nfolds = nfolds,
                            threshold = TRUE)
            fit_R <- fitVAR(X[(t + 1):n, , drop = FALSE], p = d_i,
                            nlambda = nlambda, nfolds = nfolds,
                            threshold = TRUE)
            sse <- sse + .truncated_sse(X[1:t, , drop = FALSE],     fit_L$A, d)
            sse <- sse + .truncated_sse(X[(t + 1):n, , drop = FALSE], fit_R$A, d)
        }
        obj_vals[k] <- sse
    }
    elapsed <- (proc.time() - start)[3]
    cp <- skip + which.min(obj_vals) - 1

    list(est.cp = cp, obj_vals = obj_vals, runtime = unname(elapsed))
}

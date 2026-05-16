### Fast lasso VAR fitter with a fixed regularization parameter.
###
### Used in place of sparsevar::fitVAR to skip cross-validation when running
### the simulation grid. The penalty is chosen from the manuscript formula
### lambda = 4 c0 sqrt((log p + log d) / tau), translated to glmnet's
### 1/(2n)-scaled objective.
###
### Returned list mimics sparsevar::fitVAR's interface: $A is a list of
### p x p transition matrices (one per lag) and $residuals is a matrix of
### one-step-ahead residuals.

suppressPackageStartupMessages(library(glmnet))

fast_fitVAR <- function(data, p = 1, c0 = 1.0, lambda = NULL, ...) {
    data <- as.matrix(data)
    n <- nrow(data); k <- ncol(data)
    if (n <= p + 1) stop("not enough rows to fit VAR(", p, ")")

    nrows <- n - p
    X1 <- matrix(0, nrows, k * p)
    for (q in seq_len(p)) {
        X1[, (q - 1) * k + 1:k] <- data[(p + 1 - q):(n - q), , drop = FALSE]
    }
    Y <- data[(p + 1):n, , drop = FALSE]

    if (is.null(lambda)) {
        ## Translate paper lambda = 4 c0 sqrt((log p + log d)/tau) to the
        ## glmnet scale (objective uses 1/(2n) prefactor).
        paper_lambda <- 4 * c0 * sqrt((log(max(k, nrows)) + log(p)) / nrows)
        lambda <- paper_lambda / (2 * nrows)
    }

    A_list <- vector("list", p)
    for (q in seq_len(p)) A_list[[q]] <- matrix(0, k, k)
    Yhat <- matrix(0, nrows, k)

    for (j in seq_len(k)) {
        fit <- glmnet::glmnet(X1, Y[, j], alpha = 1, lambda = lambda,
                              standardize = FALSE, intercept = FALSE)
        co <- as.numeric(coef(fit))[-1]      # drop intercept (forced to 0)
        for (q in seq_len(p)) {
            A_list[[q]][j, ] <- co[(q - 1) * k + 1:k]
        }
        Yhat[, j] <- X1 %*% co
    }
    residuals <- Y - Yhat
    list(A = A_list, residuals = residuals, lambda = lambda)
}

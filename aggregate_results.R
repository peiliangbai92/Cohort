### Aggregate simulation results into a single summary tibble.
### Produces the numbers needed to fill Table 4 (S*/M*) and Table 5 (G*) of
### the manuscript.

repo_dir <- "/Users/peiliangbai/Documents/GitHub/Cohort"
setwd(repo_dir)

scenario_meta <- list(
    S1 = list(T = 200, p = 20, M = 10, taus = 100),
    S2 = list(T = 300, p = 80, M = 10, taus = 150),
    S3 = list(T = 200, p = 20, M = 30, taus = 100),
    S4 = list(T = 200, p = 20, M = 10, taus = 100),
    S5 = list(T = 200, p = 20, M = 10, taus = 100),
    S6 = list(T = 200, p = 20, M = 10, taus = 100),
    S7 = list(T = 200, p = 20, M = 20, taus = 100),
    M1 = list(T = 300, p = 20, M = 10, taus = c(100, 200)),
    M2 = list(T = 300, p = 80, M = 10, taus = c(100, 200)),
    M3 = list(T = 1000, p = 20, M = 10, taus = c(200, 400, 600, 800)),
    M4 = list(T = 300, p = 20, M = 25, taus = c(100, 200)),
    M5 = list(T = 300, p = 20, M = 10, taus = c(100, 200))
)

select_window <- function(taus, T) {
    bounds <- numeric(length(taus) * 2)
    aug <- c(0, taus, T)
    for (j in seq_along(taus)) {
        left  <- taus[j] - (aug[j+1] - aug[j]) / 5
        right <- taus[j] + (aug[j+2] - aug[j+1]) / 5
        bounds[2*j - 1] <- left
        bounds[2*j]     <- right
    }
    bounds
}

summarise <- function(name, meta) {
    rdata <- if (name == "S1") "simulations/S1/estimated_cp.RData" else
             if (grepl("^S", name)) sprintf("simulations/%s/estimated_cp.RData", name) else
             sprintf("simulations/%s/estimated_cps.RData", name)
    if (!file.exists(rdata)) {
        return(data.frame(scenario = name, cp = NA, truth = NA,
                          mean = NA, sd = NA, pct = NA))
    }
    env <- new.env()
    load(rdata, envir = env)
    cp_obj <- if (exists("cp_result", envir = env)) get("cp_result", envir = env) else
              if (exists("est_cp", envir = env)) get("est_cp", envir = env) else NULL
    if (is.null(cp_obj)) return(NULL)

    taus <- meta$taus
    T <- meta$T
    bounds <- select_window(taus, T)

    if (grepl("^S", name)) {
        ## single CP: cp_obj is a numeric vector of length niter
        if (is.list(cp_obj)) cp_obj <- unlist(cp_obj)
        cps <- cp_obj
        successes <- cps >= bounds[1] & cps <= bounds[2]
        valid <- cps[successes]
        data.frame(scenario = name, cp = 1, truth = taus[1] / T,
                   mean = mean(valid) / T,
                   sd   = sd(valid) / T,
                   pct  = 100 * mean(successes))
    } else {
        ## multi-CP: cp_obj is a list of vectors (one per epoch)
        m0 <- length(taus)
        out <- data.frame()
        for (j in 1:m0) {
            lj <- bounds[2*j - 1]
            uj <- bounds[2*j]
            captured <- numeric()
            n_hit <- 0
            for (e in seq_along(cp_obj)) {
                cps_e <- cp_obj[[e]]
                in_band <- cps_e[cps_e >= lj & cps_e <= uj]
                if (length(in_band) >= 1) {
                    n_hit <- n_hit + 1
                    ## pick the closest one to truth if multiple hit the band
                    captured <- c(captured, in_band[which.min(abs(in_band - taus[j]))])
                }
            }
            out <- rbind(out, data.frame(
                scenario = name, cp = j, truth = taus[j] / T,
                mean = if (length(captured)) mean(captured) / T else NA,
                sd   = if (length(captured)) sd(captured) / T   else NA,
                pct  = 100 * n_hit / length(cp_obj)
            ))
        }
        out
    }
}

results <- do.call(rbind, lapply(names(scenario_meta),
                                  function(n) summarise(n, scenario_meta[[n]])))

results[, c("truth","mean","sd")] <- lapply(results[, c("truth","mean","sd")],
                                            function(x) round(x, 4))
results$pct <- round(results$pct, 1)
print(results, row.names = FALSE)
saveRDS(results, file = "logs/summary_table.rds")
write.csv(results, file = "logs/summary_table.csv", row.names = FALSE)
cat("\nSaved to logs/summary_table.{rds,csv}\n")

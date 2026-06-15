### Run all three comparison methods (SBS, BSS, fvar.seg) in parallel.
suppressPackageStartupMessages({ library(parallel) })

repo <- getwd()   # run from the repo root
setwd(repo)
log_dir <- file.path(repo, "logs")
dir.create(log_dir, showWarnings = FALSE)

scripts <- c("compare_SBS.R", "compare_BSS.R", "compare_fvarseg.R")

run_one <- function(script) {
    cmd <- sprintf("cd %s && Rscript %s > %s 2>&1",
                   shQuote(file.path(repo, "comparison")),
                   shQuote(script),
                   shQuote(file.path(log_dir, paste0(tools::file_path_sans_ext(script), ".log"))))
    t0 <- Sys.time()
    status <- system(cmd, intern = FALSE)
    dur <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    list(script = script, status = status, minutes = dur)
}

cat(sprintf("Launching %d comparison scripts at %s\n",
            length(scripts), format(Sys.time(), "%H:%M:%S")))
results <- mclapply(scripts, run_one, mc.cores = length(scripts))
cat("\n=== Summary ===\n")
for (r in results) cat(sprintf("  %s  status=%s  %.1f min\n", r$script, r$status, r$minutes))
saveRDS(results, file.path(log_dir, "compare_summary.rds"))
cat("Done at", format(Sys.time(), "%H:%M:%S"), "\n")

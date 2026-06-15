### Master parallel runner: rerun every S/M/G simulation scenario after the
### bug fixes (off-by-one in single.cp.detect, bounds.2 index, etc.).
###
### Each scenario is launched as its own R subprocess via mclapply so they run
### in parallel up to `n_cores`. Each subprocess just sources the canonical
### scenario script and writes its RData outputs into the scenario's folder.

suppressPackageStartupMessages({
    library(parallel)
})

repo_dir <- getwd()   # run from the repo root
setwd(repo_dir)

scenarios <- list(
    list(name = "S1", script = "S1/S1.R"),
    list(name = "S2", script = "S2/S2.R"),
    list(name = "S3", script = "S3/S3.R"),
    list(name = "S4", script = "S4/S4.R"),
    list(name = "S5", script = "S5/S5.R"),
    list(name = "S6", script = "S6/S6.R"),
    list(name = "S7", script = "S7/S7.R"),
    list(name = "M1", script = "M1/M1.R"),
    list(name = "M2", script = "M2/M2.R"),
    list(name = "M3", script = "M3/M3.R"),
    list(name = "M4", script = "M4/M4.R"),
    list(name = "M5", script = "M5/M5.R"),
    list(name = "G1", script = "G1/G1.R"),
    list(name = "G2", script = "G2/G2.R"),
    list(name = "G3", script = "G3/G3.R"),
    list(name = "G4", script = "G4/G4.R")
)

log_dir <- file.path(repo_dir, "logs")
dir.create(log_dir, showWarnings = FALSE)

run_one <- function(scn) {
    script_path <- file.path("simulations", scn$script)
    work_dir <- dirname(file.path(repo_dir, "simulations", scn$script))
    log_file <- file.path(log_dir, sprintf("%s.log", scn$name))
    cmd <- sprintf(
        "cd %s && Rscript %s > %s 2>&1",
        shQuote(work_dir),
        shQuote(basename(scn$script)),
        shQuote(log_file)
    )
    t0 <- Sys.time()
    status <- system(cmd, intern = FALSE)
    dur <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    list(name = scn$name, status = status, minutes = dur,
         log = log_file)
}

n_cores <- min(8, length(scenarios))
cat(sprintf("Launching %d scenarios on %d cores at %s\n",
            length(scenarios), n_cores, format(Sys.time(), "%H:%M:%S")))

results <- mclapply(scenarios, run_one, mc.cores = n_cores)

cat("\n=== Run summary ===\n")
for (r in results) {
    cat(sprintf("  %s  status=%s  %.1f min  log=%s\n",
                r$name, r$status, r$minutes, r$log))
}
saveRDS(results, file = file.path(log_dir, "run_all_summary.rds"))
cat("\nDone at", format(Sys.time(), "%H:%M:%S"), "\n")

### Application Section 5: end-to-end driver.
###
### Usage:
###   Rscript application/run_application.R          # uses cached intermediates
###   Rscript application/run_application.R --force  # rebuilds everything
###
### Pipeline (cached files in application/cache/):
###   01_preprocess_bdf.R  : raw_data/*.bdf       -> preproc_S??.rds
###   02_build_intervals.R : preproc_S??.rds      -> intervals.rds
###   03_cohort_cp.R       : intervals.rds        -> sse_lists.rds + Table 8, Fig 1
###   04_ica_clustering.R  : sse_lists.rds        -> ica.rds + Tables 9, 10, 11
###   05_networks.R        : intervals + sse      -> per_subject_A.rds + Figs 2-4
###   06_hamming.R         : per_subject_A.rds    -> hamming.rds + Tables 12-17, Fig 5
###   07_age_boxplot.R     : ica + ages CSV       -> Fig 6, permutation t-test
###   08_outlier_plot.R    : intervals + sse      -> Fig 7

args  <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
skip_age <- "--skip-age" %in% args

source("application/01_preprocess_bdf.R")
source("application/02_build_intervals.R")
source("application/03_cohort_cp.R")
source("application/04_ica_clustering.R")
source("application/05_networks.R")
source("application/06_hamming.R")
source("application/07_age_boxplot.R")
source("application/08_outlier_plot.R")

cat("\n############### Section 5.1: preprocess + intervals ###############\n")
intervals <- build_intervals(force = force)

cat("\n############### Section 5.2: cohort change-point detection ##########\n")
fits <- run_cohort_cp(force = force)
write_table8(fits)
write_fig1(fits)

cat("\n############### Section 5.2 / 5.3: ICA clustering ###################\n")
ica <- run_ica(fits)
write_table9 (ica)
write_table10(ica, fits)
write_table11(ica)
saveRDS(ica, "application/cache/ica.rds")

cat("\n############### Section 5.2: per-cluster networks ###################\n")
per_subj <- fit_per_subject_networks(intervals, fits, force = force)
write_network_figs(per_subj, ica)

cat("\n############### Section 5.2 / 5.3: Hamming distances ################\n")
H <- compute_hamming_tables(per_subj, ica)
write_hamming_tables(H)
write_fig5(H, setting = "mild")
saveRDS(H, "application/cache/hamming.rds")

cat("\n############### Section 5.3: age boxplot (Fig 6) ####################\n")
if (skip_age) {
    cat("[skip] --skip-age supplied; not running age analysis.\n")
} else {
    tryCatch(write_fig6(ica = ica),
             error = function(e) {
                 cat("[fig6 SKIPPED]:", conditionMessage(e), "\n")
                 cat("Re-run with the ages CSV filled in.\n")
             })
}

cat("\n############### Section 5.4: outlier vs typical (Fig 7) #############\n")
write_fig7(intervals = intervals, fits = fits)

cat("\nDone. Tables and figures in application/output/\n")

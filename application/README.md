# Application Section (Section 5) — reproduction

End-to-end driver for the EEG application of the cohort change-point paper.

## Inputs

* `raw_data/EEG_Cat_Study4_Resting_S{1..22}.bdf` — the 22 raw BDF recordings
  from Trujillo (2017), `doi:10.18738/T8/EG0LJI`. Already present in the repo
  (`.gitignored`).
* `application/data/subject_ages.csv` (optional) — two columns `subject, age`
  for the Section 5.3 age boxplot (Fig 6). A stub is written automatically on
  first run if missing; fill it in from the rsed2017 participant metadata.

## One-shot

```bash
Rscript application/run_application.R
```

Outputs land in `application/output/`. Intermediate RDS caches live in
`application/cache/`; pass `--force` to rebuild from scratch.

If you don't have the ages CSV yet, run with `--skip-age` to skip Fig 6.

## Individual stages (each is also runnable standalone)

| Step | Script                            | Produces                                          |
| ---- | --------------------------------- | ------------------------------------------------- |
| 1    | `01_preprocess_bdf.R`             | per-subject 21-channel detrended alpha-band RDS   |
| 2    | `02_build_intervals.R`            | 3 × 22 windowed intervals (T=1000)                |
| 3    | `03_cohort_cp.R`                  | per-subject SSEs, cohort CPs → Table 8 + Fig 1    |
| 4    | `04_ica_clustering.R`             | mild + conservative ICA → Tables 9, 10, 11        |
| 5    | `05_networks.R`                   | per-subject sparse VAR(1) → Figs 2, 3, 4          |
| 6    | `06_hamming.R`                    | Hamming tables 12–17 + Fig 5 heatmaps             |
| 7    | `07_age_boxplot.R`                | Fig 6 + permutation t-test                        |
| 8    | `08_outlier_plot.R`               | Fig 7 (subject 13 outlier vs typical)             |

## Method choices (where the paper underspecifies)

* **Channels** — paper says "21 EEG channels (NAS removed)" out of 72. We use
  the standard 10-20 montage (19 EEG + 2 mastoids = 21) and never include NAS.
  See `TARGET_CHANNELS` in `01_preprocess_bdf.R`.
* **Detrending** — iterative weighted polynomial of order 3 (de Cheveigné
  2018 style), 4 iterations, MAD-based outlier reweighting.
* **Alpha bandpass** — 4th-order Butterworth at 8-13 Hz, zero-phase via
  `signal::filtfilt`.
* **Intervals** — ±4000 samples around each empirical CP at the original
  256 Hz, subsampled by 8 → T = 1000, with the empirical CP at sample 500.
* **CP estimator** — for each subject, sparse VAR(1) fit via FISTA on both
  sides of each candidate t (`L()` in `auxiliary_functions.R`); per-subject
  CP = argmin of SSE; cohort CP = argmin of `sum_i SSE_i`.
* **ICA threshold** — mild `gamma = 3 log(p) log(n)`; conservative
  `gamma = 1.5 log(p) log(n)`. These reproduce the paper's "two clusterings"
  qualitatively; tune in `04_ica_clustering.R` if needed.
* **Network/Hamming threshold** — top 20% of `|A_ij|` across all per-subject
  matrices (one threshold per interval).

## Runtime

Roughly 40-50 minutes end-to-end on a single core. Most of the time is the
22 × 3 = 66 calls to `L()` in step 3 (each ~30-45 s).

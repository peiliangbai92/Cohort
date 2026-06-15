# Application (Section 5) — EEG reproduction

End-to-end driver for the EEG application of the cohort change-point paper.
`run_application.R` is the **paper-montage pipeline** that produced the figures
and tables in the manuscript. Run everything from the **repository root**.

## Inputs

* `application/cache/intervals.rds` — the **processed EEG** (3 × 22 montage
  windows, ~10 MB) that the method consumes. **Ships in the repo**, so the
  application reproduces without the raw recordings.
* `raw_data/EEG_Cat_Study4_Resting_S{1..22}.bdf` — the 22 raw BDF recordings
  from Trujillo (2017), `doi:10.18738/T8/EG0LJI`. Gitignored; needed **only** if
  you want to rebuild `intervals.rds` from scratch (delete the shipped one first).
* `application/data/subject_ages.csv` (optional) — two columns `subject, age`
  for the age boxplot (Fig 6). A stub with `NA` ages ships in the repo; fill it
  in from the rsed2017 participant metadata to reproduce Fig 6.

## One-shot

```bash
Rscript application/run_application.R   # figures + tables from intervals.rds
Rscript application/07_age_boxplot.R    # Fig 6 (needs subject_ages.csv filled)
```

`run_application.R` uses the shipped `intervals.rds` (skipping the raw-BDF
preprocessing); it recomputes the cohort change points and per-subject networks
each run.

`run_application.R` produces, on the paper montage:

* **Tables** → `application/output/`
  * `table8.csv` — cohort change-point summary (per-subject mean/std + cohort est.)
  * `table9.csv`, `table10.csv`, `table11.csv` — mild/conservative ICA clusterings
  * `hamming_{mild,conservative}_cp{1,2,3}.csv` — the six Hamming tables (12–17)
* **Figures** → `tex/Cohort_cp_detection_manuscript/images/`
  * `hist_cps.pdf` (Fig 1), `cp{1,2,3}_cluster*_{left,right}.pdf` (Figs 2–4 networks),
    `cp{1,2,3}_{left,right}_hamming.pdf` (Fig 5 heatmaps),
    `subject_13.pdf` / `subject_6.pdf` (Fig 7 outlier vs typical)
* **Caches** → `application/cache/` (gitignored). Stages are skipped when their
  cache exists; delete `application/cache/` (or pass `force = TRUE` to the stage
  functions) to rebuild from scratch.

`07_age_boxplot.R` writes `application/output/fig6_age_boxplot.pdf`; the
manuscript copy is `images/age_group_boxplots.pdf` (copy/rename it in).

## Pipeline stages

`run_application.R` sources stages 01–06 for their helper functions and runs the
montage pipeline on top of them. Each stage is also runnable standalone, but
**standalone 01–06 use the original montage (19 EEG + 2 mastoids); the paper's
canonical results come from `run_application.R`'s montage.**

| Stage | Script                  | Role                                                   |
| ----- | ----------------------- | ------------------------------------------------------ |
| 1     | `01_preprocess_bdf.R`   | per-subject detrended alpha-band RDS                   |
| 2     | `02_build_intervals.R`  | 3 × 22 windowed intervals (T = 1000)                   |
| 3     | `03_cohort_cp.R`        | per-subject SSEs, cohort CPs (`write_table8`)          |
| 4     | `04_ica_clustering.R`   | mild + conservative ICA (`write_table9/10/11`)         |
| 5     | `05_networks.R`         | per-subject sparse VAR(1) network estimation           |
| 6     | `06_hamming.R`          | Hamming tables + heatmaps                              |
| 7     | `07_age_boxplot.R`      | Fig 6 + permutation t-test (standalone)                |

## Method choices (where the paper underspecifies)

* **Channels (paper montage)** — 19 standard 10-20 electrodes plus the midline
  FPZ and OZ, **no mastoids** (21 total). See `TARGET_CHANNELS_V2` in
  `run_application.R`. (`TARGET_CHANNELS` in `01_preprocess_bdf.R` is the older
  19 + M1/M2 montage used only by the standalone stages.)
* **Detrending** — iterative weighted polynomial of order 3 (de Cheveigné 2018
  style), 4 iterations, MAD-based outlier reweighting.
* **Alpha bandpass** — 4th-order Butterworth at 8-13 Hz, zero-phase via
  `signal::filtfilt`.
* **Intervals** — ±4000 samples around each empirical CP at the original 256 Hz,
  subsampled by 8 → T = 1000, with the empirical CP at sample 500.
* **CP estimator** — per subject, sparse VAR(1) fit via FISTA on both sides of
  each candidate t; per-subject CP = argmin SSE; cohort CP = argmin Σᵢ SSEᵢ.
* **ICA threshold** — mild `gamma = 200`, conservative `gamma = 375` (see the
  driver's main block; tune if needed).
* **Network/Hamming threshold** — top 20% of `|A_ij|` across all per-subject
  matrices (one threshold per interval).

## Runtime

Roughly 40-50 minutes end-to-end on a single core the first time; most of it is
the 22 × 3 sparse-VAR fits. Subsequent runs reuse `application/cache/`.

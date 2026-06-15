# Cohort change-point detection

Code and data accompanying the paper on **cohort change-point detection** — a
clustering-based framework for detecting change points that are shared across a
cohort of related samples (here, multivariate time series fit with sparse
VAR models). This repository reproduces every simulation, comparison, and
the EEG application reported in the paper.

> All commands below are run **from the repository root**. The R scripts derive
> their working directory from `getwd()`, so launching elsewhere will not find
> the relative `source()` / data paths.

## Requirements

* R ≥ 4.0
* CRAN packages:
  `MTS`, `VARDetect`, `sparsevar`, `glmnet`, `mvtnorm`, `hdbinseg`,
  `signal`, `edfReader`, `igraph`, `pheatmap`, `cluster`, `factoextra`,
  `gridExtra` (plus base `parallel`, `grid`).
* `fvarseg` — factor-adjusted VAR segmentation baseline; install from its source
  repository if not on CRAN. Only needed for the `fvar.seg` comparison.

```r
install.packages(c("MTS","VARDetect","sparsevar","glmnet","mvtnorm","hdbinseg",
                   "signal","edfReader","igraph","pheatmap","cluster",
                   "factoextra","gridExtra"))
```

## Repository layout

| Path | Purpose |
| ---- | ------- |
| `auxiliary_functions.R`   | Core estimators (sparse VAR via `fista.sparse`, cohort CP, ICA clustering, Hamming). |
| `simulation_script.R`     | Data generators + simulation helpers (overrides `fitVAR` with `fast_var.R`). |
| `fast_var.R`              | Fast lasso VAR fitter used by the simulations. |
| `run_all.R`               | Runs all simulation scenarios in parallel. |
| `run_comparison.R`        | Runs the SBS / BSS / fvar.seg comparison baselines. |
| `aggregate_results.R`, `aggregate_G.R`, `aggregate_comparison.R` | Turn raw run outputs into the paper's summary tables. |
| `regenerate_table4.R`     | Emits the Table 4 LaTeX rows from the aggregated summary. |
| `simulations/`            | One folder per scenario (`S1`–`S7`, `M1`–`M5`, `G1`–`G4`); each holds its script and finished `*.RData` results. |
| `comparison/`             | The three baseline scripts and their `*.RData` results. |
| `application/`            | End-to-end EEG pipeline (see `application/README.md`). The processed EEG ships as `application/cache/intervals.rds`. |
| `data/`                   | Processed per-subject EEG segments as CSVs (supplementary; the pipeline itself runs from `intervals.rds`). |
| `tex/Cohort_cp_detection_manuscript/` | `main.tex`, `ref.bib`, and the figures in `images/`. |
| `raw_data/`               | The 22 raw EEG `.bdf` recordings (gitignored, **optional** — only needed to rebuild from scratch; see §3). |
| `logs/`                   | Run logs and aggregated summaries (gitignored, regenerable). |

## 1. Simulations → Tables 4 & 5

```bash
Rscript run_all.R            # S1–S7, M1–M5, G1–G4 → simulations/<scn>/*.RData + logs/<scn>.log
Rscript aggregate_results.R  # → logs/summary_table.{rds,csv}   (Table 4: S*/M*)
Rscript aggregate_G.R        # → logs/summary_G.{rds,csv}       (Table 5: G*)
Rscript regenerate_table4.R  # → Table 4 LaTeX rows on stdout
```

The finished `simulations/<scn>/*.RData` outputs are committed, so you can run
the three aggregate steps directly to regenerate the tables **without** the
(multi-hour) `run_all.R` re-run.

## 2. Comparison baselines → Tables 7 & 8

```bash
Rscript run_comparison.R        # SBS, BSS, fvar.seg → comparison/compare_*.RData + logs/compare_*.log
Rscript aggregate_comparison.R  # → logs/comparison_summary.{rds,csv}
```

The `compare_*.RData` results are committed; `aggregate_comparison.R` runs on
whatever is present in `comparison/`.

## 3. EEG application → Figures 1–7 and Tables 8–17

The processed EEG (`application/cache/intervals.rds`, ~10 MB) ships in the repo,
so the application reproduces **without downloading the raw recordings**:

```bash
Rscript application/run_application.R   # paper-montage pipeline (~40 min)
Rscript application/07_age_boxplot.R    # Fig 6 (needs subject ages, see below)
```

`run_application.R` reads `intervals.rds`, recomputes the cohort change points
and per-subject networks, and writes the publication figures to
`tex/Cohort_cp_detection_manuscript/images/` and the tables to
`application/output/` (Table 8 cohort CP, Tables 9–11 ICA, Tables 12–17 Hamming).

For **Fig 6**, fill in `application/data/subject_ages.csv` (`subject, age`) from
the study's participant metadata, then run `07_age_boxplot.R` (it writes
`application/output/fig6_age_boxplot.pdf`; the manuscript copy is
`images/age_group_boxplots.pdf`).

**Rebuilding from the raw BDFs (optional).** To redo preprocessing from scratch,
download the 22 recordings (Trujillo 2017, `doi:10.18738/T8/EG0LJI`) into
`raw_data/` as `EEG_Cat_Study4_Resting_S{1..22}.bdf`, delete
`application/cache/intervals.rds`, and rerun. See `application/README.md` for
per-step detail and the method choices made where the paper underspecifies.

## Manuscript

`tex/Cohort_cp_detection_manuscript/main.tex` (figures in `images/`,
bibliography in `ref.bib`).

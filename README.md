# cellrangercheck

Parse and QC **Cell Ranger** `count` / `multi` outputs across many GEM wells.

Point it at a set of run directories — each a separate GEM reaction — and get one
tidy, comparable table of pipeline QC metrics, threshold-based flagging, and a
side-by-side HTML report. Handles whichever library types are present
(Gene Expression, Antibody Capture, VDJ-T, …).

## Install

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("shihanli92/cellrangercheck")

# or with pak
# pak::pak("shihanli92/cellrangercheck")
```

Or from a local clone:

```r
# from the package directory
devtools::install("cellrangercheck")
```

Core parsing needs `readr`, `dplyr`, `tidyr`, `tibble`, `purrr`, `rlang`,
`hdf5r` and `Matrix`. The report/plots use the suggested packages `rmarkdown`,
`knitr`, `ggplot2` and `DT`.

## Quick start

```r
library(cellrangercheck)

runs <- c("path/to/sample0", "path/to/sample1")   # each a run dir or its outs/

# 1. One tidy long table across every well and sample
long <- parse_cellranger_runs(runs)

# 2. Comparable per-sample table (one row per GEM well/sample)
metrics_wide(long)

# 3. Flag metrics against thresholds and roll up per run
flags  <- flag_qc(long)                 # or flag_qc(long, my_thresholds)
qc_status(flags)                        # pass / warn / fail per run

# 4. Metadata (Cell Ranger version, chemistry, #libraries)
run_metadata(runs[1])

# 5. h5-derived metrics
barcode_ranks(runs[1])                  # barcode-rank curve + knee/inflection
umi_gene_distribution(runs)             # per-cell UMI & gene summaries
read_cr_matrix(runs[1], feature_type = "Gene Expression")  # sparse matrix

# 6. Full side-by-side HTML report (timed progress printed as it builds)
qc_report(runs, "qc_report.html")
# qc_report(runs, "qc_report.html", progress = TRUE)   # force progress off-interactive
```

`qc_report()` prints timed, per-run progress while it works, e.g.:

```
[   0.0s] Parsing 2 runs
[   0.0s]   (1/2) sample0: metrics + metadata
[   0.3s] Flagging metrics against thresholds
[   0.3s]   (1/2) sample0: barcode-rank knee
[   5.5s] Per-cell UMI / gene distributions
[   6.2s] Rendering HTML
[   7.4s] Done -> qc_report.html
```

Progress defaults to on in interactive sessions; set `progress = FALSE` to
silence it or `progress = TRUE` to force it (e.g. in scripts/logs).

## Customising thresholds

`default_thresholds()` returns an editable tibble keyed by `library_type` and
`metric` (with `category`/`grouped_by` to target the right row). Retune it and
pass it to `flag_qc()` / `qc_report()`:

```r
th <- default_thresholds()
th$fail[th$metric == "Median genes per cell"] <- 800
flag_qc(long, th)
```

## What it reads

- `outs/per_sample_outs/<sample>/metrics_summary.csv` (long, `multi`) or
  `outs/metrics_summary.csv` (wide, standalone `count`).
- `outs/config.csv` and h5 root attributes for run metadata.
- `*_feature_bc_matrix.h5` (filtered per sample, raw per GEM well) via `hdf5r`.
```

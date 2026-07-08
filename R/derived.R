#' Expected multiplet rate from recovered cell count
#'
#' Uses the 10x Genomics rule of thumb of roughly 0.8% multiplets per 1,000
#' cells recovered.
#'
#' @param n_cells Number of recovered cells.
#' @param rate_per_1k Percent multiplets per 1,000 cells. Defaults to 0.8.
#' @return Expected multiplet rate, in percent.
#' @export
expected_multiplet_rate <- function(n_cells, rate_per_1k = 0.8) {
  rate_per_1k * n_cells / 1000
}

#' Derived QC metrics as tidy long rows
#'
#' Computes the h5-derived QC metrics (mitochondrial/ribosomal content, library
#' complexity, cell-quality pass-rate, ambient RNA fraction, expected multiplet
#' rate) and returns them in the same tidy long schema as [parse_cellranger()],
#' tagged `category = "Derived"`. Bind onto the parsed metrics with
#' `dplyr::bind_rows()` and everything downstream — [metrics_wide()],
#' [flag_qc()], [assay_metric_data()] — picks them up automatically.
#'
#' @param dirs Character vector of run directories.
#' @param ids Optional run identifiers, same length as `dirs`.
#' @param ... Passed to [cell_qc()] (e.g. `min_genes`, `max_pct_mito`).
#'
#' @return A long tibble with columns matching [parse_cellranger()]:
#'   `run_id, sample_id, pipeline, category, library_type, grouped_by,
#'   group_name, metric, value, unit, value_raw`.
#' @export
derived_metrics <- function(dirs, ids = NULL, ...) {
  cq <- cell_qc(dirs, ids = ids, ...)
  amb <- ambient_fraction(dirs, ids = ids)

  # Per-sample metrics.
  per_sample <- tibble::tibble(
    run_id = rep(cq$run_id, 6),
    sample_id = rep(cq$sample_id, 6),
    metric = rep(c("Median % mitochondrial", "% cells over mito threshold",
                   "Median % ribosomal", "Median library complexity",
                   "% cells passing filter", "Expected multiplet rate"),
                 each = nrow(cq)),
    value = c(cq$median_pct_mito, cq$pct_cells_high_mito, cq$median_pct_ribo,
              cq$median_complexity, cq$pct_cells_pass,
              expected_multiplet_rate(cq$n_cells)),
    unit = rep(c("percent", "percent", "percent", "number",
                 "percent", "percent"), each = nrow(cq))
  )

  # Per-well (ambient) metric, attached to each sample of that well.
  amb_join <- dplyr::left_join(
    cq[c("run_id", "sample_id")], amb[c("run_id", "pct_ambient")],
    by = "run_id"
  )
  per_well <- tibble::tibble(
    run_id = amb_join$run_id,
    sample_id = amb_join$sample_id,
    metric = "% ambient RNA",
    value = amb_join$pct_ambient,
    unit = "percent"
  )

  out <- dplyr::bind_rows(per_sample, per_well)
  tibble::tibble(
    run_id = out$run_id,
    sample_id = out$sample_id,
    pipeline = NA_character_,
    category = "Derived",
    library_type = "Gene Expression",
    grouped_by = NA_character_,
    group_name = NA_character_,
    metric = out$metric,
    value = out$value,
    unit = out$unit,
    value_raw = as.character(round(out$value, 4))
  )
}

#' Pivot the tidy metrics table to one comparable row per sample
#'
#' Builds the side-by-side comparison table: one row per (run, sample) with a
#' column per `library_type` x `metric`. By default only the per-sample cell
#' metrics (`Category == "Cells"`, plus standalone-count rows which have no
#' category) are used, since those describe the biology being compared across
#' GEM wells. Library-level sequencing metrics remain available in the long
#' table.
#'
#' @param long A tidy tibble from [parse_cellranger()] /
#'   [parse_cellranger_runs()].
#' @param category Which `category` rows to include. Rows with a missing
#'   category (standalone count runs) are always included.
#'
#' @return A wide tibble keyed by `run_id`, `sample_id`, `pipeline` with numeric
#'   metric columns named `<library_type>__<metric>`.
#' @export
#' @examples
#' \dontrun{
#' metrics_wide(parse_cellranger_runs(c("sample0", "sample1")))
#' }
metrics_wide <- function(long, category = c("Cells", "Derived")) {
  keep <- long$category %in% category | is.na(long$category)
  sub <- long[keep, , drop = FALSE]
  # Library-level rows are keyed by group; for the per-sample comparison we
  # only keep ungrouped rows (grouped_by is NA for the "Cells" category).
  sub <- sub[is.na(sub$grouped_by), , drop = FALSE]

  sub$col <- paste0(snakecase(sub$library_type), "__", snakecase(sub$metric))
  sub <- sub[!duplicated(sub[c("run_id", "sample_id", "col")]), , drop = FALSE]

  tidyr::pivot_wider(
    sub[c("run_id", "sample_id", "pipeline", "col", "value")],
    names_from = "col",
    values_from = "value"
  )
}

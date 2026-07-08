#' Parse a wide-format (standalone count) metrics_summary.csv
#'
#' A standalone `cellranger count` run writes a single header row of metric
#' names and a single data row of values. These are mapped onto the same tidy
#' long schema used for multi runs, with `library_type = "Gene Expression"`.
#'
#' @param path Path to a wide `metrics_summary.csv`.
#' @param sample_id Sample identifier to attach.
#' @return A tidy long tibble (minus `run_id`/`pipeline`, added by the caller).
#' @keywords internal
#' @noRd
parse_metrics_count <- function(path, sample_id) {
  raw <- readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    progress = FALSE
  )
  if (nrow(raw) != 1) {
    stop("Expected a single data row in wide metrics_summary.csv '", path, "'.",
         call. = FALSE)
  }
  metric <- names(raw)
  value_raw <- unlist(raw[1, ], use.names = FALSE)
  coerced <- coerce_metric_value(value_raw)
  tibble::tibble(
    sample_id = sample_id,
    category = NA_character_,
    library_type = "Gene Expression",
    grouped_by = NA_character_,
    group_name = NA_character_,
    metric = metric,
    value = coerced$value,
    unit = coerced$unit,
    value_raw = value_raw
  )
}

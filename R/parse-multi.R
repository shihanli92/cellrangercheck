#' Parse a long-format (multi) metrics_summary.csv
#'
#' A `cellranger multi` run writes a long metrics file with columns
#' `Category, Library Type, Grouped By, Group Name, Metric Name, Metric Value`.
#' The same metric name recurs under different contexts, so the unique key is
#' all five dimensions.
#'
#' @param path Path to a per-sample `metrics_summary.csv`.
#' @param sample_id Sample identifier to attach.
#' @return A tidy long tibble (see [parse_cellranger()] for columns), minus the
#'   `run_id`/`pipeline` columns which the caller adds.
#' @keywords internal
#' @noRd
parse_metrics_multi <- function(path, sample_id) {
  raw <- readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    progress = FALSE
  )
  expected <- c("Category", "Library Type", "Grouped By", "Group Name",
                "Metric Name", "Metric Value")
  if (!all(expected %in% names(raw))) {
    stop("Unexpected multi metrics_summary.csv columns in '", path, "'.",
         call. = FALSE)
  }
  coerced <- coerce_metric_value(raw[["Metric Value"]])
  tibble::tibble(
    sample_id = sample_id,
    category = raw[["Category"]],
    library_type = raw[["Library Type"]],
    grouped_by = dplyr::na_if(raw[["Grouped By"]], ""),
    group_name = dplyr::na_if(raw[["Group Name"]], ""),
    metric = raw[["Metric Name"]],
    value = coerced$value,
    unit = coerced$unit,
    value_raw = raw[["Metric Value"]]
  )
}

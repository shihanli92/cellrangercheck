#' Flag cross-well metric outliers
#'
#' Across GEM wells, flags values that deviate from the cohort by more than
#' `nmads` median absolute deviations (a robust z-score). Useful for spotting
#' the one well that is off on cells recovered, saturation, mito content, etc.
#'
#' @param wide A wide metrics table from [metrics_wide()].
#' @param nmads Threshold in MADs. Defaults to 3.
#' @param metrics Optional character vector of metric columns to consider;
#'   defaults to all numeric columns.
#'
#' @return A long tibble of the flagged values: `run_id`, `sample_id`, `metric`,
#'   `value`, `median`, `mad`, `n_mads`, `is_outlier`. Columns with fewer than 3
#'   non-missing wells or zero MAD are skipped (no robust call possible). Returns
#'   zero rows when nothing is flagged.
#' @export
flag_outliers <- function(wide, nmads = 3, metrics = NULL) {
  key_cols <- intersect(c("run_id", "sample_id", "pipeline"), names(wide))
  num_cols <- names(wide)[vapply(wide, is.numeric, logical(1))]
  if (!is.null(metrics)) num_cols <- intersect(num_cols, metrics)

  out <- purrr::map_dfr(num_cols, function(col) {
    v <- wide[[col]]
    ok <- !is.na(v)
    if (sum(ok) < 3) return(NULL)
    med <- stats::median(v, na.rm = TRUE)
    md <- stats::mad(v, na.rm = TRUE)
    if (is.na(md) || md == 0) return(NULL)
    n_mads <- abs(v - med) / md
    tibble::tibble(
      run_id = wide$run_id,
      sample_id = if ("sample_id" %in% key_cols) wide$sample_id else NA_character_,
      metric = col,
      value = v,
      median = med,
      mad = md,
      n_mads = n_mads,
      is_outlier = !is.na(n_mads) & n_mads > nmads
    )
  })
  if (is.null(out) || nrow(out) == 0) {
    return(tibble::tibble(
      run_id = character(0), sample_id = character(0), metric = character(0),
      value = numeric(0), median = numeric(0), mad = numeric(0),
      n_mads = numeric(0), is_outlier = logical(0)
    ))
  }
  out[out$is_outlier %in% TRUE, , drop = FALSE]
}

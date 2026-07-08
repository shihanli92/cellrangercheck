#' Flag metrics against QC thresholds
#'
#' Joins the tidy metrics table to a threshold table and assigns each matched
#' metric a status of `"pass"`, `"warn"` or `"fail"`.
#'
#' @param long A tidy metrics tibble from [parse_cellranger()] /
#'   [parse_cellranger_runs()].
#' @param thresholds A threshold tibble; see [default_thresholds()].
#'
#' @return The matched rows with added columns `direction`, `warn`, `fail` and
#'   `status`. Use [qc_status()] to roll up to one status per run/sample.
#' @export
#' @examples
#' \dontrun{
#' flags <- flag_qc(parse_cellranger_runs(c("sample0", "sample1")))
#' qc_status(flags)
#' }
flag_qc <- function(long, thresholds = default_thresholds()) {
  # Match on library_type + metric, then filter to the intended category /
  # grouping so an ambiguous metric name (e.g. "Cells") hits the right row.
  joined <- dplyr::inner_join(
    long, thresholds,
    by = c("library_type", "metric"),
    relationship = "many-to-many"
  )
  match_ok <-
    (is.na(joined$category.y) | joined$category.x == joined$category.y) &
    (is.na(joined$grouped_by.y) | joined$grouped_by.x == joined$grouped_by.y)
  joined <- joined[match_ok, , drop = FALSE]

  joined$status <- qc_status_one(joined$value, joined$direction,
                                 joined$warn, joined$fail)

  out <- joined[c("run_id", "sample_id", "library_type", "metric",
                  "category.x", "grouped_by.x", "group_name",
                  "value", "unit", "direction", "warn", "fail", "status")]
  names(out)[names(out) == "category.x"] <- "category"
  names(out)[names(out) == "grouped_by.x"] <- "grouped_by"
  tibble::as_tibble(out)
}

#' Vectorised single-metric status
#' @keywords internal
#' @noRd
qc_status_one <- function(value, direction, warn, fail) {
  status <- rep(NA_character_, length(value))
  is_min <- direction == "min"
  is_max <- direction == "max"

  # Higher-is-better
  status[is_min & value >= warn] <- "pass"
  status[is_min & value < warn & value >= fail] <- "warn"
  status[is_min & value < fail] <- "fail"

  # Lower-is-better
  status[is_max & value <= warn] <- "pass"
  status[is_max & value > warn & value <= fail] <- "warn"
  status[is_max & value > fail] <- "fail"

  status[is.na(value)] <- NA_character_
  status
}

#' Roll flagged metrics up to one status per run and sample
#'
#' @param flags Output of [flag_qc()].
#' @return A tibble with `run_id`, `sample_id`, counts of `n_pass`, `n_warn`,
#'   `n_fail`, and the worst-case `qc_status`.
#' @export
qc_status <- function(flags) {
  worst <- function(s) {
    s <- s[!is.na(s)]
    if (length(s) == 0) return(NA_character_)
    if ("fail" %in% s) return("fail")
    if ("warn" %in% s) return("warn")
    "pass"
  }
  dplyr::summarise(
    dplyr::group_by(flags, .data$run_id, .data$sample_id),
    n_pass = sum(.data$status == "pass", na.rm = TRUE),
    n_warn = sum(.data$status == "warn", na.rm = TRUE),
    n_fail = sum(.data$status == "fail", na.rm = TRUE),
    qc_status = worst(.data$status),
    .groups = "drop"
  )
}

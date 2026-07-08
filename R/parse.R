#' Parse a single Cell Ranger run into a tidy metrics table
#'
#' Detects whether `dir` is a `cellranger multi` or `cellranger count` run and
#' returns every reported metric in one long tibble. Multi runs contribute one
#' block of rows per sample under `per_sample_outs/`.
#'
#' @param dir Path to a Cell Ranger run directory (the run root or its `outs/`).
#' @param id Optional run identifier. Defaults to the run directory's basename.
#'
#' @return A tibble with columns `run_id`, `sample_id`, `pipeline`, `category`,
#'   `library_type`, `grouped_by`, `group_name`, `metric`, `value`, `unit` and
#'   `value_raw`.
#' @seealso [parse_cellranger_runs()], [metrics_wide()]
#' @export
#' @examples
#' \dontrun{
#' parse_cellranger("sample0")
#' }
parse_cellranger <- function(dir, id = NULL) {
  info <- detect_pipeline(dir)
  run_id <- id %||% run_id_from_dir(dir)

  rows <- purrr::pmap_dfr(
    list(info$samples$metrics_csv, info$samples$sample_id),
    function(csv, sid) {
      if (info$pipeline == "multi") {
        parse_metrics_multi(csv, sid)
      } else {
        parse_metrics_count(csv, sid)
      }
    }
  )

  tibble::add_column(
    rows,
    run_id = run_id,
    pipeline = info$pipeline,
    .before = "sample_id"
  )
}

#' Parse many Cell Ranger runs into one combined tidy table
#'
#' The headline entry point: point it at a vector of run directories (each a
#' separate GEM reaction) and get a single long tibble spanning every well and
#' sample, ready for [metrics_wide()] and [flag_qc()].
#'
#' @param dirs Character vector of run directories.
#' @param ids Optional character vector of run identifiers, the same length as
#'   `dirs`. Defaults to each directory's basename.
#'
#' @return A combined long tibble; see [parse_cellranger()] for columns.
#' @export
#' @examples
#' \dontrun{
#' long <- parse_cellranger_runs(c("sample0", "sample1"))
#' }
parse_cellranger_runs <- function(dirs, ids = NULL) {
  if (!is.null(ids) && length(ids) != length(dirs)) {
    stop("`ids` must be the same length as `dirs`.", call. = FALSE)
  }
  ids <- ids %||% rep(NA_character_, length(dirs))
  purrr::pmap_dfr(list(dirs, ids), function(d, i) {
    parse_cellranger(d, id = if (is.na(i)) NULL else i)
  })
}

#' Derive a run id from a run directory path
#' @keywords internal
#' @noRd
run_id_from_dir <- function(dir) {
  dir <- normalizePath(dir, mustWork = TRUE)
  if (basename(dir) == "outs") {
    return(basename(dirname(dir)))
  }
  basename(dir)
}

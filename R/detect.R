#' Resolve the `outs/` directory for a Cell Ranger run
#'
#' Accepts either a run root (containing `outs/`) or the `outs/` directory
#' itself.
#'
#' @param dir Path to a run directory.
#' @return Absolute path to the `outs/` directory.
#' @keywords internal
#' @noRd
resolve_outs <- function(dir) {
  dir <- normalizePath(dir, mustWork = TRUE)
  if (basename(dir) == "outs" && dir.exists(dir)) {
    return(dir)
  }
  outs <- file.path(dir, "outs")
  if (dir.exists(outs)) {
    return(normalizePath(outs))
  }
  # Perhaps the user pointed directly at a per-sample dir or the outs itself
  # lacks the name; fall back to the dir if it looks like an outs.
  if (dir.exists(file.path(dir, "per_sample_outs")) ||
      file.exists(file.path(dir, "metrics_summary.csv"))) {
    return(dir)
  }
  stop("Could not find an 'outs/' directory under '", dir, "'.", call. = FALSE)
}

#' Detect the Cell Ranger pipeline and locate metric files
#'
#' Classifies a run as `multi` (long-format metrics_summary.csv, one per sample
#' under `per_sample_outs/`) or `count` (a single wide metrics_summary.csv).
#'
#' @param dir Path to a run directory (root or its `outs/`).
#' @return A list with `pipeline` (`"multi"` or `"count"`), `outs` (path), and
#'   `samples`: a tibble with one row per sample carrying `sample_id`,
#'   `metrics_csv` and `count_dir` (the per-sample `count/` folder, or the run
#'   `outs/` for a standalone count run).
#' @keywords internal
#' @noRd
detect_pipeline <- function(dir) {
  outs <- resolve_outs(dir)
  per_sample <- file.path(outs, "per_sample_outs")

  if (dir.exists(per_sample)) {
    sample_dirs <- list.dirs(per_sample, recursive = FALSE)
    samples <- purrr::map_dfr(sample_dirs, function(sd) {
      tibble::tibble(
        sample_id = basename(sd),
        metrics_csv = file.path(sd, "metrics_summary.csv"),
        count_dir = file.path(sd, "count")
      )
    })
    samples <- samples[file.exists(samples$metrics_csv), , drop = FALSE]
    if (nrow(samples) == 0) {
      stop("Found 'per_sample_outs/' but no metrics_summary.csv under '",
           outs, "'.", call. = FALSE)
    }
    return(list(pipeline = "multi", outs = outs, samples = samples))
  }

  # Standalone count run: a single wide metrics_summary.csv in outs/.
  metrics_csv <- file.path(outs, "metrics_summary.csv")
  if (file.exists(metrics_csv)) {
    return(list(
      pipeline = "count",
      outs = outs,
      samples = tibble::tibble(
        sample_id = basename(dirname(outs)),
        metrics_csv = metrics_csv,
        count_dir = outs
      )
    ))
  }

  stop("No metrics_summary.csv found for run '", dir, "'.", call. = FALSE)
}

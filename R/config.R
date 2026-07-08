#' Read run-level metadata for a Cell Ranger run
#'
#' Combines the multi config (`outs/config.csv`, when present) with attributes
#' read from a feature-barcode HDF5 matrix (Cell Ranger version and chemistry).
#'
#' @param dir Path to a run directory (root or its `outs/`).
#'
#' @return A one-row tibble with `run_id`, `software_version`,
#'   `chemistry_description`, and `n_libraries` (from `config.csv`, or `NA`).
#' @export
run_metadata <- function(dir) {
  info <- detect_pipeline(dir)
  run_id <- run_id_from_dir(dir)

  # Config: number of libraries declared in the [libraries] section.
  cfg_path <- file.path(info$outs, "config.csv")
  n_libraries <- NA_integer_
  if (file.exists(cfg_path)) {
    n_libraries <- count_config_libraries(cfg_path)
  }

  # h5 attributes from the first available matrix.
  h5_path <- find_matrix_h5(info, which = "filtered")
  attrs <- list(software_version = NA_character_,
                chemistry_description = NA_character_)
  if (!is.na(h5_path)) {
    attrs <- read_h5_run_attrs(h5_path)
  }

  tibble::tibble(
    run_id = run_id,
    software_version = attrs$software_version,
    chemistry_description = attrs$chemistry_description,
    n_libraries = n_libraries
  )
}

#' Count declared libraries in a multi config.csv
#' @keywords internal
#' @noRd
count_config_libraries <- function(cfg_path) {
  lines <- readLines(cfg_path, warn = FALSE)
  lines <- trimws(lines)
  sect <- NA_character_
  n <- 0L
  header_seen <- FALSE
  for (ln in lines) {
    if (grepl("^\\[.*\\]$", ln)) {
      sect <- gsub("^\\[|\\]$", "", ln)
      header_seen <- FALSE
      next
    }
    if (identical(sect, "libraries") && nzchar(ln)) {
      if (!header_seen) {
        header_seen <- TRUE  # first line is the column header
      } else {
        n <- n + 1L
      }
    }
  }
  if (n == 0L) NA_integer_ else n
}

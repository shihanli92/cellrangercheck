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

#' Parse the [libraries] section of a multi config.csv
#'
#' @param cfg_path Path to a run's `outs/config.csv`.
#' @return A tibble with columns `fastq_id`, `fastqs`, `feature_types` (one row
#'   per declared library), or `NULL` if the section is absent/empty.
#' @keywords internal
#' @noRd
parse_config_libraries <- function(cfg_path) {
  lines <- readLines(cfg_path, warn = FALSE)
  lines <- trimws(lines)
  sect <- NA_character_
  header <- NULL
  rows <- list()
  for (ln in lines) {
    if (grepl("^\\[.*\\]$", ln)) {
      sect <- gsub("^\\[|\\]$", "", ln)
      next
    }
    if (identical(sect, "libraries") && nzchar(ln)) {
      fields <- trimws(strsplit(ln, ",", fixed = TRUE)[[1]])
      if (is.null(header)) {
        header <- fields          # first line is the column header
      } else {
        rows[[length(rows) + 1L]] <- fields
      }
    }
  }
  if (is.null(header) || length(rows) == 0L) {
    return(NULL)
  }
  mat <- do.call(rbind, lapply(rows, function(r) {
    length(r) <- length(header)   # pad short rows
    r
  }))
  df <- tibble::as_tibble(as.data.frame(mat, stringsAsFactors = FALSE))
  names(df) <- header
  # Normalise the columns we rely on downstream.
  std <- c(fastq_id = "fastq_id", fastqs = "fastqs",
           feature_types = "feature_types")
  for (nm in names(std)) {
    if (!nm %in% names(df)) df[[nm]] <- NA_character_
  }
  df[c("fastq_id", "fastqs", "feature_types")]
}

#' Count declared libraries in a multi config.csv
#' @keywords internal
#' @noRd
count_config_libraries <- function(cfg_path) {
  libs <- parse_config_libraries(cfg_path)
  if (is.null(libs)) NA_integer_ else nrow(libs)
}

#' Build a fastq_id -> reaction map from run config.csv files
#'
#' Reads the `[libraries]` section of each run's `outs/config.csv` and returns
#' the fastq prefixes declared for that run, so FastQC outputs (named by fastq
#' prefix) can be attributed to the GEM reaction that produced them.
#'
#' @param dirs Character vector of run directories.
#' @param ids Optional run identifiers, same length as `dirs`.
#'
#' @return A tibble with `run_id`, `fastq_id`, `fastqs`, `feature_types`. Runs
#'   without a readable `config.csv` contribute no rows (with a warning).
#' @export
library_map <- function(dirs, ids = NULL) {
  ids <- ids %||% rep(NA_character_, length(dirs))
  empty <- tibble::tibble(
    run_id = character(0), fastq_id = character(0),
    fastqs = character(0), feature_types = character(0)
  )
  out <- purrr::pmap_dfr(list(dirs, ids), function(d, i) {
    info <- detect_pipeline(d)
    run_id <- if (is.na(i)) run_id_from_dir(d) else i
    cfg_path <- file.path(info$outs, "config.csv")
    if (!file.exists(cfg_path)) {
      warning("No config.csv for run '", run_id, "'; cannot map its fastqs.",
              call. = FALSE)
      return(NULL)
    }
    libs <- parse_config_libraries(cfg_path)
    if (is.null(libs)) return(NULL)
    tibble::add_column(libs, run_id = run_id, .before = "fastq_id")
  })
  if (nrow(out) == 0) empty else out
}

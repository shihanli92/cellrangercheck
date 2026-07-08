#' Render a side-by-side QC report across GEM wells
#'
#' Parses each run, flags metrics, and renders an HTML report comparing all GEM
#' reactions. Requires the suggested packages `rmarkdown`, `knitr`, `ggplot2`
#' and `DT`.
#'
#' The metrics (including the slower h5 reads) are computed up front so progress
#' can be reported; the Rmd template only renders the precomputed results.
#'
#' @param dirs Character vector of run directories.
#' @param output_file Path for the rendered HTML file.
#' @param ids Optional run identifiers, same length as `dirs`.
#' @param thresholds Threshold tibble; see [default_thresholds()].
#' @param include_h5 Whether to compute h5-derived metrics (barcode-rank knee,
#'   UMI/gene distributions). Set `FALSE` to skip if matrices are absent.
#' @param fastqc Optional character vector of directories to scan for FastQC
#'   outputs. When supplied, a FastQC section (per-module PASS/WARN/FAIL grid and
#'   basic statistics, grouped by reaction) is appended. FastQC files are matched
#'   to reactions via each run's `config.csv` fastq prefixes.
#' @param progress Whether to print timed progress messages while building the
#'   report. Defaults to `TRUE` in interactive sessions.
#'
#' @return The path to the rendered report, invisibly.
#' @export
qc_report <- function(dirs, output_file, ids = NULL,
                      thresholds = default_thresholds(), include_h5 = TRUE,
                      fastqc = NULL, progress = interactive()) {
  rlang::check_installed(c("rmarkdown", "knitr", "ggplot2", "DT"))
  template <- system.file("rmd", "report.Rmd", package = "cellrangercheck")
  if (!nzchar(template)) {
    stop("Report template not found in the installed package.", call. = FALSE)
  }
  # Make the output path absolute even if the file does not exist yet:
  # rmarkdown::render interprets a relative output_file relative to the
  # template's directory, so a bare filename would be written next to the
  # template instead of the caller's working directory.
  output_file <- file.path(
    normalizePath(dirname(output_file), mustWork = FALSE),
    basename(output_file)
  )
  # Rendering runs in the template's directory, so resolve run paths to
  # absolute first or relative `dirs` will not be found.
  dirs <- normalizePath(dirs, mustWork = TRUE)
  run_ids <- ids %||% basename(dirs)
  n <- length(dirs)

  step <- progress_stepper(progress)

  # ---- Parse metrics + metadata, per run so progress is visible ----
  step(sprintf("Parsing %d run%s", n, if (n == 1) "" else "s"))
  long_list <- vector("list", n)
  meta_list <- vector("list", n)
  for (i in seq_len(n)) {
    step(sprintf("  (%d/%d) %s: metrics + metadata", i, n, run_ids[i]))
    long_list[[i]] <- parse_cellranger(dirs[i], id = run_ids[i])
    m <- run_metadata(dirs[i])
    m$run_id <- run_ids[i]
    meta_list[[i]] <- m
  }
  long <- dplyr::bind_rows(long_list)
  meta <- dplyr::bind_rows(meta_list)

  step("Flagging metrics against thresholds")
  wide <- metrics_wide(long)
  flags <- flag_qc(long, thresholds)
  status <- qc_status(flags)
  status$overall <- status$qc_status  # updated below if FastQC is included

  # ---- h5-derived metrics (the slow part) ----
  br_df <- NULL
  dist <- NULL
  if (isTRUE(include_h5)) {
    step("Reading h5 matrices (barcode ranks)")
    br_list <- vector("list", n)
    for (i in seq_len(n)) {
      step(sprintf("  (%d/%d) %s: barcode-rank knee", i, n, run_ids[i]))
      br <- tryCatch(barcode_ranks(dirs[i]), error = function(e) NULL)
      if (!is.null(br)) {
        d <- br$ranks
        d$run_id <- run_ids[i]
        d$knee <- br$knee
        br_list[[i]] <- d
      }
    }
    br_df <- dplyr::bind_rows(br_list)
    if (nrow(br_df) == 0) br_df <- NULL

    step("Per-cell UMI / gene distributions")
    dist <- tryCatch(umi_gene_distribution(dirs, ids = run_ids),
                     error = function(e) NULL)
  }

  # ---- optional FastQC section ----
  fqc <- NULL
  if (!is.null(fastqc)) {
    step("Scanning + parsing FastQC outputs")
    fqc <- tryCatch(fastqc_by_reaction(fastqc, dirs, ids = run_ids),
                    error = function(e) {
                      warning(conditionMessage(e), call. = FALSE); NULL
                    })
    if (is.null(fqc)) {
      step("  no FastQC outputs found")
    } else {
      # Fold the FastQC roll-up into the overall per-reaction status.
      fqs <- fastqc_status(fqc)
      status <- dplyr::left_join(status, fqs, by = "run_id")
      status$overall <- mapply(worst_status, status$qc_status,
                               status$fastqc_status, USE.NAMES = FALSE)
    }
  }

  step("Rendering HTML")
  rmarkdown::render(
    template,
    output_file = output_file,
    params = list(
      precomputed = list(
        long = long, meta = meta, wide = wide,
        flags = flags, status = status, br = br_df, dist = dist,
        fqc = fqc
      ),
      include_h5 = include_h5
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
  step(sprintf("Done -> %s", output_file))
  invisible(output_file)
}

#' Build a timed progress-message emitter
#'
#' Returns a function that prints each message prefixed with the elapsed time
#' since the stepper was created. A no-op when `enabled` is `FALSE`.
#'
#' @param enabled Logical; whether to print.
#' @return A function taking a single message string.
#' @keywords internal
#' @noRd
progress_stepper <- function(enabled = TRUE) {
  start <- Sys.time()
  function(msg) {
    if (!enabled) return(invisible())
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    message(sprintf("[%6.1fs] %s", elapsed, msg))
    invisible()
  }
}

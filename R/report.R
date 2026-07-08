#' Render a side-by-side QC report across GEM wells
#'
#' Parses each run, flags metrics, and renders an HTML report comparing all GEM
#' reactions. Requires the suggested packages `rmarkdown`, `knitr`, `ggplot2`
#' and `DT`.
#'
#' @param dirs Character vector of run directories.
#' @param output_file Path for the rendered HTML file.
#' @param ids Optional run identifiers, same length as `dirs`.
#' @param thresholds Threshold tibble; see [default_thresholds()].
#' @param include_h5 Whether to compute h5-derived metrics (barcode-rank knee,
#'   UMI/gene distributions). Set `FALSE` to skip if matrices are absent.
#'
#' @return The path to the rendered report, invisibly.
#' @export
qc_report <- function(dirs, output_file, ids = NULL,
                      thresholds = default_thresholds(), include_h5 = TRUE) {
  rlang::check_installed(c("rmarkdown", "knitr", "ggplot2", "DT"))
  template <- system.file("rmd", "report.Rmd", package = "cellrangercheck")
  if (!nzchar(template)) {
    stop("Report template not found in the installed package.", call. = FALSE)
  }
  output_file <- normalizePath(output_file, mustWork = FALSE)
  # Rendering runs in the template's directory, so resolve run paths to
  # absolute first or relative `dirs` will not be found.
  dirs <- normalizePath(dirs, mustWork = TRUE)
  rmarkdown::render(
    template,
    output_file = output_file,
    params = list(
      dirs = as.list(dirs),
      ids = if (is.null(ids)) NULL else as.list(ids),
      thresholds = thresholds,
      include_h5 = include_h5
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
  invisible(output_file)
}

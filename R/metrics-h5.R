#' Barcode-rank curve and knee/inflection estimates for a run
#'
#' Reads the raw (unfiltered) feature-barcode matrix, ranks barcodes by total
#' UMI count, and estimates the knee and inflection points on the log-log
#' barcode-rank curve. This is a lightweight sanity check on cell calling; it
#' does not reproduce the exact `cellranger`/`emptyDrops` algorithm.
#'
#' @param dir Path to a run directory (root or its `outs/`).
#' @param feature_type Feature type used for the UMI totals. Defaults to
#'   `"Gene Expression"`.
#'
#' @return A list with `ranks` (a tibble of `rank`, `total_umi` for the
#'   deduplicated curve), `knee` and `inflection` (estimated UMI thresholds),
#'   and `n_barcodes_above_knee`.
#' @export
barcode_ranks <- function(dir, feature_type = "Gene Expression") {
  mat <- read_cr_matrix(dir, feature_type = feature_type, which = "raw")
  totals <- Matrix::colSums(mat)
  totals <- sort(totals[totals > 0], decreasing = TRUE)

  ranks <- tibble::tibble(rank = seq_along(totals), total_umi = as.numeric(totals))
  # Deduplicate on total_umi (using the last rank at each level) for a stable
  # curve, mirroring the standard barcode-rank presentation.
  dedup <- ranks[!duplicated(ranks$total_umi, fromLast = TRUE), , drop = FALSE]

  kp <- knee_point(dedup$rank, dedup$total_umi)

  list(
    ranks = dedup,
    knee = kp$knee,
    inflection = kp$inflection,
    n_barcodes_above_knee = sum(totals >= kp$knee)
  )
}

#' Estimate knee and inflection on a log-log barcode-rank curve
#'
#' The inflection is taken as the point of steepest descent (most negative
#' first derivative) of log10(total) vs log10(rank); the knee as the point of
#' maximum curvature (most negative second derivative). Both are returned as UMI
#' totals.
#'
#' @param rank Integer ranks (ascending).
#' @param total Total UMI counts (descending), same length as `rank`.
#' @return A list with `knee` and `inflection` UMI totals.
#' @keywords internal
#' @noRd
knee_point <- function(rank, total) {
  if (length(rank) < 4) {
    return(list(knee = NA_real_, inflection = NA_real_))
  }
  x <- log10(rank)
  y <- log10(total)
  d1 <- diff(y) / diff(x)
  inflection_idx <- which.min(d1) + 1L
  d2 <- diff(d1) / diff(x[-1])
  knee_idx <- which.min(d2) + 1L
  list(
    knee = total[min(knee_idx, length(total))],
    inflection = total[min(inflection_idx, length(total))]
  )
}

#' Per-cell UMI and gene distributions across runs
#'
#' Reads each run's filtered (cell) matrix and summarises per-cell total UMIs
#' and detected genes.
#'
#' @param dirs Character vector of run directories.
#' @param ids Optional run identifiers, same length as `dirs`.
#' @param feature_type Feature type to summarise. Defaults to
#'   `"Gene Expression"`.
#'
#' @return A tibble with one row per (run, sample): `run_id`, `sample_id`,
#'   `n_cells`, `median_umi`, `mad_umi`, `q25_umi`, `q75_umi`,
#'   `median_genes`, `mad_genes`, `q25_genes`, `q75_genes`.
#' @export
umi_gene_distribution <- function(dirs, ids = NULL, feature_type = "Gene Expression") {
  ids <- ids %||% rep(NA_character_, length(dirs))
  purrr::pmap_dfr(list(dirs, ids), function(d, i) {
    info <- detect_pipeline(d)
    run_id <- if (is.na(i)) run_id_from_dir(d) else i
    purrr::map_dfr(info$samples$sample_id, function(sid) {
      parsed <- read_cr_matrix_features(d, feature_type = feature_type,
                                        sample_id = sid, which = "filtered")
      pc <- per_cell_qc(parsed$mat, parsed$features)
      tibble::tibble(
        run_id = run_id,
        sample_id = sid,
        n_cells = nrow(pc),
        median_umi = stats::median(pc$umi),
        mad_umi = stats::mad(pc$umi),
        q25_umi = stats::quantile(pc$umi, 0.25, names = FALSE),
        q75_umi = stats::quantile(pc$umi, 0.75, names = FALSE),
        median_genes = stats::median(pc$genes),
        mad_genes = stats::mad(pc$genes),
        q25_genes = stats::quantile(pc$genes, 0.25, names = FALSE),
        q75_genes = stats::quantile(pc$genes, 0.75, names = FALSE)
      )
    })
  })
}

#' Ambient (background) RNA fraction per GEM well
#'
#' Estimates how much of the sequenced signal comes from cell-free ("empty")
#' droplets: reads the raw (unfiltered) matrix and the set of called-cell
#' barcodes, then reports the fraction of total UMIs that fall outside cells. A
#' high value flags ambient/soup contamination.
#'
#' @param dirs Character vector of run directories.
#' @param ids Optional run identifiers, same length as `dirs`.
#' @param feature_type Feature type to use. Defaults to `"Gene Expression"`.
#'
#' @return A tibble with one row per run: `run_id`, `total_umi` (raw),
#'   `cell_umi`, `pct_ambient`, `n_empty`, `mean_umi_per_empty`.
#' @export
ambient_fraction <- function(dirs, ids = NULL, feature_type = "Gene Expression") {
  ids <- ids %||% rep(NA_character_, length(dirs))
  purrr::pmap_dfr(list(dirs, ids), function(d, i) {
    info <- detect_pipeline(d)
    run_id <- if (is.na(i)) run_id_from_dir(d) else i

    raw <- read_cr_matrix(d, feature_type = feature_type, which = "raw")
    raw_totals <- Matrix::colSums(raw)
    total_umi <- sum(raw_totals)

    # Cells = union of the run's per-sample filtered barcodes (raw is well-level).
    cell_bcs <- unique(unlist(lapply(info$samples$sample_id, function(sid) {
      colnames(read_cr_matrix(d, feature_type = feature_type, sample_id = sid,
                              which = "filtered"))
    })))
    is_cell <- colnames(raw) %in% cell_bcs
    cell_umi <- sum(raw_totals[is_cell])
    n_empty <- sum(!is_cell)

    tibble::tibble(
      run_id = run_id,
      total_umi = total_umi,
      cell_umi = cell_umi,
      pct_ambient = if (total_umi > 0) 100 * (total_umi - cell_umi) / total_umi
                    else NA_real_,
      n_empty = n_empty,
      mean_umi_per_empty = if (n_empty > 0) sum(raw_totals[!is_cell]) / n_empty
                           else NA_real_
    )
  })
}

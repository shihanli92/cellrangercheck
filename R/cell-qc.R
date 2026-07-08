#' Per-cell QC metrics from a gene-expression matrix
#'
#' @param mat A genes x cells sparse matrix.
#' @param features A feature tibble with a `name` column (gene symbols).
#' @param mito_pattern,ribo_pattern Regexes (case-insensitive) matching
#'   mitochondrial / ribosomal gene symbols.
#' @return A tibble with one row per cell: `umi`, `genes`, `pct_mito`,
#'   `pct_ribo`, `complexity`.
#' @keywords internal
#' @noRd
per_cell_qc <- function(mat, features,
                        mito_pattern = "^mt-", ribo_pattern = "^rp[sl]") {
  umi <- Matrix::colSums(mat)
  genes <- Matrix::colSums(mat > 0)

  mito_rows <- grepl(mito_pattern, features$name, ignore.case = TRUE)
  ribo_rows <- grepl(ribo_pattern, features$name, ignore.case = TRUE)
  mito_umi <- if (any(mito_rows)) Matrix::colSums(mat[mito_rows, , drop = FALSE])
              else rep(0, ncol(mat))
  ribo_umi <- if (any(ribo_rows)) Matrix::colSums(mat[ribo_rows, , drop = FALSE])
              else rep(0, ncol(mat))

  safe_pct <- function(num, den) ifelse(den > 0, 100 * num / den, NA_real_)
  # log10(genes)/log10(umi): novelty score, undefined for umi<=1.
  complexity <- ifelse(umi > 1 & genes > 0, log10(genes) / log10(umi), NA_real_)

  tibble::tibble(
    umi = as.numeric(umi),
    genes = as.numeric(genes),
    pct_mito = safe_pct(mito_umi, umi),
    pct_ribo = safe_pct(ribo_umi, umi),
    complexity = complexity
  )
}

#' Per-cell QC summary per GEM well
#'
#' Computes mitochondrial and ribosomal content, library complexity and a
#' cell-quality pass-rate from each run's filtered (cell) matrix. Mitochondrial
#' and ribosomal genes are matched on gene symbol, case-insensitively, so the
#' defaults work for both mouse (`mt-`, `Rps`/`Rpl`) and human (`MT-`,
#' `RPS`/`RPL`).
#'
#' @param dirs Character vector of run directories.
#' @param ids Optional run identifiers, same length as `dirs`.
#' @param min_genes,max_pct_mito Filters defining a "usable" cell: at least
#'   `min_genes` detected genes and at most `max_pct_mito` percent mitochondrial.
#' @param feature_type Feature type to use. Defaults to `"Gene Expression"`.
#' @param mito_pattern,ribo_pattern Gene-symbol regexes (case-insensitive).
#'
#' @return A tibble with one row per (run, sample): `run_id`, `sample_id`,
#'   `n_cells`, `median_umi`, `median_genes`, `median_pct_mito`, `mad_pct_mito`,
#'   `pct_cells_high_mito`, `median_pct_ribo`, `median_complexity`,
#'   `n_cells_pass`, `pct_cells_pass`.
#' @export
cell_qc <- function(dirs, ids = NULL, min_genes = 200, max_pct_mito = 20,
                    feature_type = "Gene Expression",
                    mito_pattern = "^mt-", ribo_pattern = "^rp[sl]") {
  ids <- ids %||% rep(NA_character_, length(dirs))
  purrr::pmap_dfr(list(dirs, ids), function(d, i) {
    info <- detect_pipeline(d)
    run_id <- if (is.na(i)) run_id_from_dir(d) else i
    purrr::map_dfr(info$samples$sample_id, function(sid) {
      parsed <- read_cr_matrix_features(d, feature_type = feature_type,
                                        sample_id = sid, which = "filtered")
      pc <- per_cell_qc(parsed$mat, parsed$features,
                        mito_pattern = mito_pattern, ribo_pattern = ribo_pattern)
      pass <- pc$genes >= min_genes &
        (is.na(pc$pct_mito) | pc$pct_mito <= max_pct_mito)
      n <- nrow(pc)
      tibble::tibble(
        run_id = run_id,
        sample_id = sid,
        n_cells = n,
        median_umi = stats::median(pc$umi),
        median_genes = stats::median(pc$genes),
        median_pct_mito = stats::median(pc$pct_mito, na.rm = TRUE),
        mad_pct_mito = stats::mad(pc$pct_mito, na.rm = TRUE),
        pct_cells_high_mito = 100 * sum(pc$pct_mito > max_pct_mito,
                                        na.rm = TRUE) / n,
        median_pct_ribo = stats::median(pc$pct_ribo, na.rm = TRUE),
        median_complexity = stats::median(pc$complexity, na.rm = TRUE),
        n_cells_pass = sum(pass),
        pct_cells_pass = 100 * sum(pass) / n
      )
    })
  })
}

multi_dir <- function() system.file("extdata", "gm", package = "cellrangercheck")

test_that("per_cell_qc computes mito/ribo/complexity against a known matrix", {
  # 3 genes x 2 cells: gene 2 is mito, gene 3 is ribo.
  m <- Matrix::Matrix(matrix(c(10, 5, 5,     # cell1: umi 20, mito 5, ribo 5
                               30, 0, 10),   # cell2: umi 40, mito 0, ribo 10
                             nrow = 3), sparse = TRUE)
  feats <- tibble::tibble(name = c("GeneA", "mt-Nd1", "Rps2"))
  pc <- cellrangercheck:::per_cell_qc(m, feats)
  expect_equal(pc$umi, c(20, 40))
  expect_equal(pc$pct_mito, c(25, 0))
  expect_equal(pc$pct_ribo, c(25, 25))
  expect_equal(pc$genes, c(3, 2))
})

test_that("cell_qc summarises per well with mito/ribo/pass-rate columns", {
  cq <- cell_qc(multi_dir(), ids = "gm")
  expect_equal(cq$run_id, "gm")
  expect_equal(cq$n_cells, 5)
  expect_true(all(c("median_pct_mito", "median_pct_ribo", "median_complexity",
                    "pct_cells_pass", "pct_cells_high_mito") %in% names(cq)))
  expect_true(cq$median_pct_mito > 0)        # fixture has mt- genes with counts
  expect_gte(cq$pct_cells_pass, 0)
  expect_lte(cq$pct_cells_pass, 100)
})

test_that("pass-rate responds to the filter thresholds", {
  lax <- cell_qc(multi_dir(), min_genes = 1, max_pct_mito = 100)
  strict <- cell_qc(multi_dir(), min_genes = 1e6, max_pct_mito = 0)
  expect_equal(lax$pct_cells_pass, 100)
  expect_equal(strict$pct_cells_pass, 0)
})

test_that("umi_gene_distribution still returns its documented columns", {
  d <- umi_gene_distribution(multi_dir(), ids = "gm")
  expect_true(all(c("n_cells", "median_umi", "mad_umi", "median_genes") %in%
                    names(d)))
  expect_equal(d$n_cells, 5)
})

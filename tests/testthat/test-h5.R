multi_dir <- function() system.file("extdata", "gm", package = "cellrangercheck")

test_that("read_cr_matrix returns a sparse features x barcodes matrix", {
  m <- read_cr_matrix(multi_dir())
  expect_s4_class(m, "dgCMatrix")
  expect_equal(ncol(m), 5)      # 5 cells in the filtered fixture
  expect_equal(nrow(m), 6)      # 4 GEX + 2 antibody features
})

test_that("feature_type subsets the matrix rows", {
  gex <- read_cr_matrix(multi_dir(), feature_type = "Gene Expression")
  expect_equal(nrow(gex), 4)
  ab <- read_cr_matrix(multi_dir(), feature_type = "Antibody Capture")
  expect_equal(nrow(ab), 2)
})

test_that("barcode_ranks estimates a knee above background", {
  br <- barcode_ranks(multi_dir())
  expect_true(is.finite(br$knee))
  expect_gte(br$n_barcodes_above_knee, 1)
  expect_true(all(diff(br$ranks$total_umi) <= 0))  # descending curve
})

test_that("umi_gene_distribution summarises per-cell counts", {
  d <- umi_gene_distribution(multi_dir(), ids = "gm")
  expect_equal(d$n_cells, 5)
  expect_true(d$median_umi > 0)
  expect_true(d$median_genes > 0)
})

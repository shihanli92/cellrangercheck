test_that("flag_outliers catches a well that deviates by >3 MADs", {
  wide <- tibble::tibble(
    run_id = paste0("r", 1:5),
    sample_id = paste0("r", 1:5),
    pipeline = "multi",
    gene_expression__median_genes_per_cell = c(2000, 2010, 1990, 2005, 50)
  )
  out <- flag_outliers(wide)
  expect_equal(nrow(out), 1)
  expect_equal(out$run_id, "r5")
  expect_true(out$is_outlier)
  expect_gt(out$n_mads, 3)
})

test_that("flag_outliers needs at least 3 wells", {
  wide <- tibble::tibble(
    run_id = c("a", "b"), sample_id = c("a", "b"), pipeline = "multi",
    m = c(1, 100)
  )
  expect_equal(nrow(flag_outliers(wide)), 0)
})

test_that("flag_outliers ignores zero-MAD (constant) columns", {
  wide <- tibble::tibble(
    run_id = paste0("r", 1:4), sample_id = paste0("r", 1:4), pipeline = "m",
    constant = rep(5, 4)
  )
  expect_equal(nrow(flag_outliers(wide)), 0)
})

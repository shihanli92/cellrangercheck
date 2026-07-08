multi_dir <- function() system.file("extdata", "gm", package = "cellrangercheck")

test_that("ambient_fraction computes a plausible background fraction", {
  a <- ambient_fraction(multi_dir(), ids = "gm")
  expect_equal(a$run_id, "gm")
  expect_equal(a$n_empty, 40)                 # raw fixture has 40 empty droplets
  expect_true(a$pct_ambient >= 0 && a$pct_ambient <= 100)
  expect_true(a$cell_umi <= a$total_umi)
  # Cells dominate the signal in the fixture, so ambient is small.
  expect_lt(a$pct_ambient, 50)
})

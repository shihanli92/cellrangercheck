test_that("metric values coerce with correct units", {
  res <- coerce_metric_value(c("86.3%", "47,538", "2.0", "412.07", ""))
  expect_equal(res$value, c(86.3, 47538, 2.0, 412.07, NA))
  expect_equal(res$unit, c("percent", "count", "number", "number", NA))
})

test_that("snakecase normalises metric labels", {
  expect_equal(snakecase("Median UMI counts per cell"),
               "median_umi_counts_per_cell")
  expect_equal(snakecase("Cells With Productive V-J Spanning Pair"),
               "cells_with_productive_v_j_spanning_pair")
})

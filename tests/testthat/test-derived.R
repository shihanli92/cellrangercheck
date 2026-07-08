multi_dir <- function() system.file("extdata", "gm", package = "cellrangercheck")

test_that("expected_multiplet_rate follows the 10x heuristic", {
  expect_equal(expected_multiplet_rate(1000), 0.8)
  expect_equal(expected_multiplet_rate(5000), 4.0)
  expect_equal(expected_multiplet_rate(0), 0)
})

test_that("derived_metrics returns the tidy long schema tagged Derived", {
  d <- derived_metrics(multi_dir(), ids = "gm")
  expect_setequal(
    names(d),
    c("run_id", "sample_id", "pipeline", "category", "library_type",
      "grouped_by", "group_name", "metric", "value", "unit", "value_raw")
  )
  expect_true(all(d$category == "Derived"))
  expect_true(all(d$library_type == "Gene Expression"))
  expect_true(all(c("Median % mitochondrial", "% cells passing filter",
                    "% ambient RNA", "Expected multiplet rate") %in% d$metric))
})

test_that("derived rows bind and flow through wide + flag_qc", {
  long <- dplyr::bind_rows(
    parse_cellranger(multi_dir(), id = "gm"),
    derived_metrics(multi_dir(), ids = "gm")
  )
  w <- metrics_wide(long)
  expect_true("gene_expression__median_mitochondrial" %in% names(w))

  flags <- flag_qc(long)
  # The mito / pass-rate / ambient thresholds should now match derived rows.
  expect_true(any(flags$metric == "Median % mitochondrial"))
  expect_true(any(flags$metric == "% cells passing filter"))
})

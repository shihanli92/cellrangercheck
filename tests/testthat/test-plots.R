multi_dir <- function() system.file("extdata", "gm", package = "cellrangercheck")

test_that("assay_metric_data returns one tidy row per library-metric", {
  md <- assay_metric_data(parse_cellranger(multi_dir(), id = "gm"))
  expect_true(all(c("run_id", "sample_id", "library_type", "library",
                    "metric_label", "value", "unit") %in% names(md)))
  expect_setequal(unique(md$library_type),
                  c("Gene Expression", "Antibody Capture"))
  expect_true(all(is.finite(md$value)))
})

test_that("metric labels disambiguate cell vs sequencing metrics", {
  md <- assay_metric_data(parse_cellranger(multi_dir()))
  expect_true(any(grepl("\\[cell\\]$", md$metric_label)))
  expect_true(any(grepl("\\[seq\\]$", md$metric_label)))
  # 'Cells' appears in both categories but stays distinct via the tag.
  cells <- md$metric_label[grepl("^Cells ", md$metric_label) &
                             md$library_type == "Gene Expression"]
  expect_true("Cells [seq]" %in% cells)
})

test_that("assay_metric_plot builds a plotly widget", {
  skip_if_not_installed("plotly")
  md <- assay_metric_data(parse_cellranger(multi_dir()))
  p <- assay_metric_plot(md, "Gene Expression")
  expect_s3_class(p, "plotly")
})

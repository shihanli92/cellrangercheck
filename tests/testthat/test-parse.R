multi_dir <- function() system.file("extdata", "gm", package = "cellrangercheck")
count_dir <- function() system.file("extdata", "gc", package = "cellrangercheck")

test_that("multi and count runs are detected", {
  expect_equal(detect_pipeline(multi_dir())$pipeline, "multi")
  expect_equal(detect_pipeline(count_dir())$pipeline, "count")
})

test_that("parse_cellranger returns the tidy long schema", {
  long <- parse_cellranger(multi_dir(), id = "gm")
  expect_setequal(
    names(long),
    c("run_id", "pipeline", "sample_id", "category", "library_type",
      "grouped_by", "group_name", "metric", "value", "unit", "value_raw")
  )
  expect_equal(unique(long$run_id), "gm")
  expect_equal(unique(long$sample_id), "S1")
  expect_true(all(c("Gene Expression", "Antibody Capture") %in% long$library_type))
})

test_that("the five context dimensions keep same-named metrics distinct", {
  long <- parse_cellranger(multi_dir())
  cells_rows <- long[long$library_type == "Gene Expression" &
                       long$metric == "Cells", ]
  # "Cells" appears under both the Cells and Library categories.
  expect_true(all(c("Cells", "Library") %in% cells_rows$category))
  expect_gt(nrow(cells_rows), 1)
})

test_that("wide count metrics parse and coerce", {
  long <- parse_cellranger(count_dir())
  expect_equal(long$pipeline[1], "count")
  est <- long$value[long$metric == "Estimated Number of Cells"]
  expect_equal(est, 1234)
})

test_that("parse_cellranger_runs binds multiple wells with correct ids", {
  long <- parse_cellranger_runs(c(multi_dir(), count_dir()),
                                ids = c("A", "B"))
  expect_setequal(unique(long$run_id), c("A", "B"))
  expect_setequal(unique(long$pipeline), c("multi", "count"))
})

test_that("metrics_wide pivots to one row per sample", {
  long <- parse_cellranger(multi_dir())
  w <- metrics_wide(long)
  expect_equal(nrow(w), 1)
  expect_true("gene_expression__median_genes_per_cell" %in% names(w))
})

test_that("run_metadata reads version and chemistry", {
  meta <- run_metadata(multi_dir())
  expect_equal(meta$software_version, "cellranger-9.0.1")
  expect_match(meta$chemistry_description, "Single Cell")
  expect_equal(meta$n_libraries, 2L)
})

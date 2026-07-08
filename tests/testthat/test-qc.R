multi_dir <- function() system.file("extdata", "gm", package = "cellrangercheck")

test_that("qc_status_one classifies min/max directions", {
  expect_equal(
    qc_status_one(c(95, 80, 60), "min", warn = 90, fail = 75),
    c("pass", "warn", "fail")
  )
  expect_equal(
    qc_status_one(c(5, 15, 30), "max", warn = 10, fail = 20),
    c("pass", "warn", "fail")
  )
  expect_true(is.na(qc_status_one(NA_real_, "min", 90, 75)))
})

test_that("flag_qc matches thresholds to the intended category", {
  long <- parse_cellranger(multi_dir())
  flags <- flag_qc(long)
  expect_true(all(flags$status %in% c("pass", "warn", "fail")))
  # The GEX 'Confidently mapped reads in cells' is a Cells-category metric.
  cmap <- flags[flags$library_type == "Gene Expression" &
                  flags$metric == "Confidently mapped reads in cells", ]
  expect_equal(nrow(cmap), 1)
  expect_equal(cmap$category, "Cells")
})

test_that("qc_status rolls up to the worst per sample", {
  long <- parse_cellranger(multi_dir())
  status <- qc_status(flag_qc(long))
  expect_equal(nrow(status), 1)
  expect_true(status$qc_status %in% c("pass", "warn", "fail"))
  expect_equal(status$n_pass + status$n_warn + status$n_fail,
               nrow(flag_qc(long)))
})

test_that("custom thresholds change the outcome", {
  long <- parse_cellranger(multi_dir())
  strict <- default_thresholds()
  strict$warn[strict$metric == "Median genes per cell"] <- 1e6
  strict$fail[strict$metric == "Median genes per cell"] <- 1e5
  flags <- flag_qc(long, strict)
  gpc <- flags[flags$metric == "Median genes per cell", ]
  expect_equal(gpc$status, "fail")
})

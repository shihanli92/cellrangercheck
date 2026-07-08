# End-to-end checks against real Cell Ranger runs. Skipped unless the runs are
# available (set CELLRANGERCHECK_RUNS to a comma-separated list of run dirs).
real_runs <- function() {
  env <- Sys.getenv("CELLRANGERCHECK_RUNS", "")
  if (nzchar(env)) {
    runs <- strsplit(env, ",", fixed = TRUE)[[1]]
    return(runs[dir.exists(runs)])
  }
  character(0)
}

test_that("real multi runs parse, pivot and flag end to end", {
  runs <- real_runs()
  skip_if(length(runs) < 1, "No real runs available")

  long <- parse_cellranger_runs(runs)
  expect_gt(nrow(long), 0)
  expect_true(all(c("run_id", "sample_id", "library_type", "metric",
                    "value") %in% names(long)))

  wide <- metrics_wide(long)
  expect_equal(nrow(wide), length(unique(paste(long$run_id, long$sample_id))))

  status <- qc_status(flag_qc(long))
  expect_true(all(status$qc_status %in% c("pass", "warn", "fail")))
})

test_that("real runs expose h5-derived metrics", {
  runs <- real_runs()
  skip_if(length(runs) < 1, "No real runs available")
  d <- umi_gene_distribution(runs)
  expect_true(all(d$n_cells > 0))
  br <- barcode_ranks(runs[1])
  expect_true(is.finite(br$knee))
})

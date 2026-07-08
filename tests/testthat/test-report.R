test_that("qc_report writes to the intended path (relative resolves to cwd)", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")
  skip_if_not_installed("DT")
  skip_if_not_installed("ggplot2")

  d <- system.file("extdata", "gm", package = "cellrangercheck")
  tmp <- tempfile("qcrep")
  dir.create(tmp)
  tmp <- normalizePath(tmp)
  old <- setwd(tmp)
  on.exit(setwd(old), add = TRUE)

  out <- qc_report(d, "report.html", include_h5 = TRUE, progress = FALSE)

  # Lands in the working directory, not next to the Rmd template.
  expect_true(file.exists("report.html"))
  expect_true(file.exists(out))
  expect_equal(normalizePath(out), normalizePath(file.path(tmp, "report.html")))

  tmpl_dir <- dirname(system.file("rmd", "report.Rmd",
                                  package = "cellrangercheck"))
  expect_false(file.exists(file.path(tmpl_dir, "report.html")))
})

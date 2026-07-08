fqc_dir <- function() system.file("extdata", "fastqc", package = "cellrangercheck")
gm_dir <- function() system.file("extdata", "gm", package = "cellrangercheck")

test_that("library_map reads fastq prefixes from config.csv", {
  m <- library_map(gm_dir(), ids = "gm")
  expect_setequal(m$fastq_id, c("GEX0", "SP0"))
  expect_equal(m$run_id, c("gm", "gm"))
  expect_true("Gene Expression" %in% m$feature_types)
})

test_that("parse_fastqc reads both extracted dirs and zips", {
  fq <- parse_fastqc(fqc_dir())
  expect_equal(nrow(fq$stats), 3)                 # 2 extracted + 1 zip
  expect_setequal(fq$stats$read, c("R1", "R2"))
  # zip entry (SP0) parsed for basic stats
  sp <- fq$stats[grepl("^SP0", fq$stats$filename), ]
  expect_equal(sp$total_sequences, 5e5)
  expect_equal(sp$pct_dup, 55.0)                  # 100 - 45.0 dedup
})

test_that("basic statistics are parsed from fastqc_data.txt", {
  fq <- parse_fastqc(fqc_dir())
  gex_r1 <- fq$stats[fq$stats$filename == "GEX0_S1_L001_R1_001.fastq.gz", ]
  expect_equal(gex_r1$total_sequences, 1e6)
  expect_equal(gex_r1$pct_gc, 48)
  expect_equal(gex_r1$pct_dup, 12.5)              # 100 - 87.5
  expect_equal(gex_r1$read_length, "28")
})

test_that("module statuses are captured", {
  fq <- parse_fastqc(fqc_dir())
  expect_true(all(fq$modules$status %in% c("PASS", "WARN", "FAIL")))
  ov <- fq$modules[fq$modules$filename == "GEX0_S1_L001_R1_001.fastq.gz" &
                     fq$modules$module == "Overrepresented sequences", ]
  expect_equal(ov$status, "FAIL")
})

test_that("fastqc_by_reaction attributes files to the right reaction", {
  fr <- fastqc_by_reaction(fqc_dir(), gm_dir(), ids = "gm")
  expect_true(all(fr$stats$run_id == "gm"))
  gex <- fr$stats[grepl("^GEX0", fr$stats$filename), ]
  expect_true(all(gex$feature_types == "Gene Expression"))
  sp <- fr$stats[grepl("^SP0", fr$stats$filename), ]
  expect_true(all(sp$feature_types == "Antibody Capture"))
  expect_false(any(is.na(fr$stats$run_id)))       # everything matched
})

test_that("infer_fastq_id resolves Illumina sample names and prefixes", {
  ids <- c("GEX0", "SP0")
  expect_equal(infer_fastq_id("GEX0_S1_L001_R1_001.fastq.gz", ids), "GEX0")
  expect_equal(infer_fastq_id("SP0_S2_L002_R2_001.fastq.gz", ids), "SP0")
  expect_true(is.na(infer_fastq_id("VDJ9_S1_L001_R1_001.fastq.gz", ids)))
})

test_that("unmatched FastQC files warn but do not error", {
  # gc run has no config.csv libraries mapping GEX0/SP0 prefixes
  fr <- NULL
  w <- capture_warnings(
    fr <- fastqc_by_reaction(fqc_dir(),
                             system.file("extdata", "gc",
                                         package = "cellrangercheck"))
  )
  expect_true(any(grepl("matched no fastq_id", w)))
  expect_true(all(is.na(fr$stats$run_id)))
})

test_that("parse_fastqc returns NULL when nothing is found", {
  empty <- tempfile("empty_scan")
  dir.create(empty)
  expect_null(parse_fastqc(empty))
})

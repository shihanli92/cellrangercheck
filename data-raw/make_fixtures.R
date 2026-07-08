# Generates the small synthetic fixtures under inst/extdata used by tests.
# Run from the package root: Rscript data-raw/make_fixtures.R
library(Matrix)

ext <- "inst/extdata"

# ---- tiny 10x-format h5 writer -------------------------------------------
write_10x_h5 <- function(path, mat, feature_id, feature_name, feature_type,
                         barcodes, software_version = "cellranger-9.0.1",
                         chemistry = "Single Cell 3' v3 (polyA)") {
  if (file.exists(path)) file.remove(path)
  mat <- as(mat, "CsparseMatrix")
  f <- hdf5r::H5File$new(path, mode = "w")
  on.exit(f$close_all(), add = TRUE)
  g <- f$create_group("matrix")
  g[["data"]] <- as.integer(mat@x)
  g[["indices"]] <- as.integer(mat@i)
  g[["indptr"]] <- as.integer(mat@p)
  g[["shape"]] <- as.integer(dim(mat))
  g[["barcodes"]] <- barcodes
  ft <- g$create_group("features")
  ft[["id"]] <- feature_id
  ft[["name"]] <- feature_name
  ft[["feature_type"]] <- feature_type
  ft[["_all_tag_keys"]] <- "genome"
  ft[["genome"]] <- rep("GRCh38", length(feature_id))
  hdf5r::h5attr(f, "software_version") <- software_version
  hdf5r::h5attr(f, "chemistry_description") <- chemistry
  hdf5r::h5attr(f, "filetype") <- "matrix"
  hdf5r::h5attr(f, "version") <- 2L
  invisible(path)
}

set.seed(1)  # deterministic

# Features: 4 gene-expression + 2 antibody.
feature_id <- c("ENSG1", "ENSG2", "ENSG3", "ENSG4", "AB1", "AB2")
feature_name <- c("GeneA", "GeneB", "GeneC", "GeneD", "HTO1", "HTO2")
feature_type <- c(rep("Gene Expression", 4), rep("Antibody Capture", 2))

# ---- filtered matrix: 5 clear cells ----
n_cells <- 5
filt_bc <- sprintf("CELL%02d-1", seq_len(n_cells))
filt <- matrix(rpois(length(feature_id) * n_cells, lambda = 40),
               nrow = length(feature_id))
filt <- Matrix(filt, sparse = TRUE)
write_10x_h5(
  file.path(ext, "gm/outs/per_sample_outs/S1/count",
            "sample_filtered_feature_bc_matrix.h5"),
  filt, feature_id, feature_name, feature_type, filt_bc
)

# ---- raw matrix: 5 cells with high counts + 40 empty-ish droplets ----
n_empty <- 40
raw_bc <- c(filt_bc, sprintf("BG%03d-1", seq_len(n_empty)))
big <- matrix(rpois(length(feature_id) * n_cells, lambda = 400),
              nrow = length(feature_id))
small <- matrix(rpois(length(feature_id) * n_empty, lambda = 1),
                nrow = length(feature_id))
raw <- Matrix(cbind(big, small), sparse = TRUE)
write_10x_h5(
  file.path(ext, "gm/outs/multi/count", "raw_feature_bc_matrix.h5"),
  raw, feature_id, feature_name, feature_type, raw_bc
)

# ---- standalone count fixture: wide metrics_summary.csv + matrices ----
dir.create(file.path(ext, "gc/outs"), recursive = TRUE, showWarnings = FALSE)
wide <- data.frame(
  check.names = FALSE,
  "Estimated Number of Cells" = "1,234",
  "Mean Reads per Cell" = "45,678",
  "Median Genes per Cell" = "2,001",
  "Sequencing Saturation" = "62.3%",
  "Reads Mapped Confidently to Transcriptome" = "81.0%",
  "Valid Barcodes" = "97.9%"
)
write.csv(wide, file.path(ext, "gc/outs/metrics_summary.csv"),
          row.names = FALSE, quote = TRUE)

# Give the count fixture GEX-only matrices so h5 helpers work there too.
gex_id <- feature_id[1:4]; gex_name <- feature_name[1:4]
gex_type <- feature_type[1:4]
write_10x_h5(
  file.path(ext, "gc/outs", "filtered_feature_bc_matrix.h5"),
  filt[1:4, ], gex_id, gex_name, gex_type, filt_bc
)
write_10x_h5(
  file.path(ext, "gc/outs", "raw_feature_bc_matrix.h5"),
  raw[1:4, ], gex_id, gex_name, gex_type, raw_bc
)

# ---- FastQC fixtures (match the gm config fastq_ids: GEX0, SP0) ----
fqc_root <- file.path(ext, "fastqc")
unlink(fqc_root, recursive = TRUE)
dir.create(fqc_root, recursive = TRUE, showWarnings = FALSE)

write_fastqc <- function(root, sample_name, dedup_pct, gc, total, len,
                         statuses, as_zip = FALSE, overrep = NULL) {
  fq <- paste0(sample_name, ".fastq.gz")
  base <- paste0(sample_name, "_fastqc")
  summary_txt <- paste(
    apply(statuses, 1, function(r) paste(r[["status"]], r[["module"]], fq,
                                         sep = "\t")),
    collapse = "\n")
  overrep_block <- c(">>Overrepresented sequences\tpass", ">>END_MODULE")
  if (!is.null(overrep)) {
    overrep_block <- c(
      ">>Overrepresented sequences\tfail",
      "#Sequence\tCount\tPercentage\tPossible Source",
      apply(overrep, 1, function(r) paste(r[["sequence"]], r[["count"]],
                                          r[["percentage"]], r[["source"]],
                                          sep = "\t")),
      ">>END_MODULE")
  }
  data_txt <- paste(c(
    "##FastQC\t0.12.1",
    ">>Basic Statistics\tpass",
    "#Measure\tValue",
    paste0("Filename\t", fq),
    "File type\tConventional base calls",
    "Encoding\tSanger / Illumina 1.9",
    paste0("Total Sequences\t", total),
    "Sequences flagged as poor quality\t0",
    paste0("Sequence length\t", len),
    paste0("%GC\t", gc),
    ">>END_MODULE",
    ">>Sequence Duplication Levels\tpass",
    paste0("#Total Deduplicated Percentage\t", dedup_pct),
    ">>END_MODULE",
    overrep_block
  ), collapse = "\n")

  if (as_zip) {
    tmp <- file.path(tempdir(), base)
    dir.create(tmp, showWarnings = FALSE)
    writeLines(summary_txt, file.path(tmp, "summary.txt"))
    writeLines(data_txt, file.path(tmp, "fastqc_data.txt"))
    zip_path <- file.path(normalizePath(root), paste0(base, ".zip"))
    old <- setwd(tempdir()); on.exit(setwd(old), add = TRUE)
    utils::zip(zip_path, base, flags = "-rq")
    unlink(tmp, recursive = TRUE)
  } else {
    d <- file.path(root, base)
    dir.create(d, showWarnings = FALSE)
    writeLines(summary_txt, file.path(d, "summary.txt"))
    writeLines(data_txt, file.path(d, "fastqc_data.txt"))
  }
}

gex_status <- data.frame(
  status = c("PASS", "PASS", "WARN", "PASS", "FAIL"),
  module = c("Basic Statistics", "Per base sequence quality",
             "Per base sequence content", "Per sequence GC content",
             "Overrepresented sequences"),
  stringsAsFactors = FALSE)
sp_status <- data.frame(
  status = c("PASS", "WARN", "PASS", "PASS", "WARN"),
  module = gex_status$module, stringsAsFactors = FALSE)

overrep_tbl <- data.frame(
  sequence = c("AAAAAAAAAAAAAAAAAAAAAAAA", "GGGGGGGGGGGGGGGGGGGGGGGG"),
  count = c(52000, 18000),
  percentage = c(5.2, 1.8),
  source = c("No Hit", "No Hit"),
  stringsAsFactors = FALSE)

# Extracted-dir FastQC for GEX0 (R1 + R2); zipped FastQC for SP0.
write_fastqc(fqc_root, "GEX0_S1_L001_R1_001", 87.5, 48, 1000000, 28, gex_status,
             overrep = overrep_tbl)
write_fastqc(fqc_root, "GEX0_S1_L001_R2_001", 82.1, 47, 1000000, 90, gex_status,
             overrep = overrep_tbl)
write_fastqc(fqc_root, "SP0_S2_L001_R1_001", 45.0, 51, 500000, 28, sp_status,
             as_zip = TRUE)

message("Fixtures written under ", normalizePath(ext))

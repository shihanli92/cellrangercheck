#' Discover FastQC outputs under one or more directories
#'
#' Recursively finds both zipped (`*_fastqc.zip`) and already-extracted
#' (`*_fastqc/` containing `fastqc_data.txt`) FastQC outputs.
#'
#' @param dirs Character vector of directories to scan.
#' @return A character vector of FastQC output paths (zip files and extracted
#'   directories), de-duplicated so a zip and its extracted copy are not both
#'   returned.
#' @keywords internal
#' @noRd
find_fastqc_outputs <- function(dirs) {
  dirs <- dirs[dir.exists(dirs)]
  zips <- character(0)
  exdirs <- character(0)
  for (d in dirs) {
    zips <- c(zips, list.files(d, pattern = "_fastqc\\.zip$",
                               recursive = TRUE, full.names = TRUE))
    data_txt <- list.files(d, pattern = "^fastqc_data\\.txt$",
                           recursive = TRUE, full.names = TRUE)
    exdirs <- c(exdirs, dirname(data_txt))
  }
  # Drop extracted dirs whose sibling .zip is also present (avoid duplicates).
  ex_keep <- if (length(exdirs)) {
    exdirs[!paste0(exdirs, ".zip") %in% zips]
  } else {
    character(0)
  }
  unique(c(zips, ex_keep))
}

#' Read a member file from a FastQC output (zip or extracted dir)
#' @keywords internal
#' @noRd
read_fastqc_member <- function(path, member) {
  if (dir.exists(path)) {
    fp <- file.path(path, member)
    if (!file.exists(fp)) return(NULL)
    return(readLines(fp, warn = FALSE))
  }
  # Zip archive: locate the member ending in the requested name.
  entries <- tryCatch(utils::unzip(path, list = TRUE)$Name,
                      error = function(e) character(0))
  hit <- entries[grepl(paste0("(^|/)", member, "$"), entries)]
  if (length(hit) == 0) return(NULL)
  con <- unz(path, hit[1])
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

#' Parse a single FastQC output into module statuses and basic statistics
#' @keywords internal
#' @noRd
read_fastqc_one <- function(path) {
  summary_lines <- read_fastqc_member(path, "summary.txt")
  data_lines <- read_fastqc_member(path, "fastqc_data.txt")
  if (is.null(summary_lines) && is.null(data_lines)) {
    warning("Could not read FastQC output '", path, "'.", call. = FALSE)
    return(NULL)
  }

  # summary.txt: <status>\t<module>\t<filename>
  modules <- NULL
  filename <- NA_character_
  if (!is.null(summary_lines)) {
    parts <- strsplit(summary_lines, "\t", fixed = TRUE)
    parts <- parts[lengths(parts) >= 3]
    modules <- tibble::tibble(
      module = vapply(parts, `[`, character(1), 2),
      status = toupper(vapply(parts, `[`, character(1), 1))
    )
    filename <- parts[[1]][3]
  }

  # fastqc_data.txt: Basic Statistics block + duplication percentage.
  stats <- list(total_sequences = NA_real_, pct_gc = NA_real_,
                read_length = NA_character_, pct_dup = NA_real_)
  if (!is.null(data_lines)) {
    kv <- function(key) {
      hit <- grep(paste0("^", key, "\t"), data_lines, value = TRUE)
      if (length(hit) == 0) return(NA_character_)
      sub(paste0("^", key, "\t"), "", hit[1])
    }
    if (is.na(filename)) filename <- kv("Filename")
    stats$total_sequences <- suppressWarnings(as.numeric(kv("Total Sequences")))
    stats$pct_gc <- suppressWarnings(as.numeric(kv("%GC")))
    stats$read_length <- kv("Sequence length")
    dedup <- kv("#Total Deduplicated Percentage")
    if (!is.na(dedup)) {
      stats$pct_dup <- round(100 - as.numeric(dedup), 2)
    }
  }

  list(path = path, filename = filename, modules = modules, stats = stats)
}

#' Infer the Illumina read tag (R1/R2/I1/I2) from a fastq filename
#' @keywords internal
#' @noRd
fastqc_read_tag <- function(filename) {
  m <- regmatches(filename, regexpr("_(R[12]|I[12])_", filename))
  if (length(m) == 0 || m == "") return(NA_character_)
  gsub("_", "", m)
}

#' Parse FastQC outputs found under directories
#'
#' Discovers and parses every FastQC output under `dirs`, returning per-module
#' PASS/WARN/FAIL statuses and basic statistics per fastq file.
#'
#' @param dirs Character vector of directories to scan for FastQC outputs.
#'
#' @return A list with two tibbles:
#'   * `stats`: one row per FastQC file with `path`, `filename`, `read`,
#'     `total_sequences`, `pct_gc`, `read_length`, `pct_dup`.
#'   * `modules`: long form `path`, `filename`, `read`, `module`, `status`.
#'   Returns `NULL` if no FastQC outputs are found.
#' @export
parse_fastqc <- function(dirs) {
  paths <- find_fastqc_outputs(dirs)
  if (length(paths) == 0) return(NULL)
  parsed <- lapply(paths, read_fastqc_one)
  parsed <- Filter(Negate(is.null), parsed)
  if (length(parsed) == 0) return(NULL)

  stats <- purrr::map_dfr(parsed, function(p) {
    tibble::tibble(
      path = p$path,
      filename = p$filename,
      read = fastqc_read_tag(p$filename),
      total_sequences = p$stats$total_sequences,
      pct_gc = p$stats$pct_gc,
      read_length = p$stats$read_length,
      pct_dup = p$stats$pct_dup
    )
  })
  modules <- purrr::map_dfr(parsed, function(p) {
    if (is.null(p$modules)) return(NULL)
    tibble::add_column(
      p$modules,
      path = p$path, filename = p$filename,
      read = fastqc_read_tag(p$filename), .before = "module"
    )
  })
  list(stats = stats, modules = modules)
}

#' Attribute a fastq filename to a declared fastq_id
#' @keywords internal
#' @noRd
infer_fastq_id <- function(filename, known_ids) {
  known_ids <- known_ids[!is.na(known_ids)]
  if (length(known_ids) == 0 || is.na(filename)) return(NA_character_)
  base <- basename(filename)
  # Illumina sample name = everything before the _S<n>_ token.
  token <- sub("_S[0-9]+_.*$", "", base)
  if (token %in% known_ids) return(token)
  cand <- known_ids[startsWith(base, known_ids)]
  if (length(cand)) cand[which.max(nchar(cand))] else NA_character_
}

#' Parse FastQC outputs and attribute them to Cell Ranger reactions
#'
#' Scans `fastqc_dirs` for FastQC outputs, reads each run's `outs/config.csv`
#' to learn its declared fastq prefixes, and joins the two so every FastQC file
#' is tagged with the GEM reaction (run) and library it belongs to.
#'
#' @param fastqc_dirs Character vector of directories to scan for FastQC output.
#' @param run_dirs Character vector of Cell Ranger run directories.
#' @param ids Optional run identifiers, same length as `run_dirs`.
#'
#' @return A list with `stats` and `modules` tibbles (as in [parse_fastqc()]),
#'   each gaining `run_id`, `fastq_id` and `feature_types` columns. Files that
#'   match no declared fastq_id keep `NA` in those columns (with a warning).
#'   Returns `NULL` if no FastQC outputs are found.
#' @export
fastqc_by_reaction <- function(fastqc_dirs, run_dirs, ids = NULL) {
  fq <- parse_fastqc(fastqc_dirs)
  if (is.null(fq)) return(NULL)
  map <- library_map(run_dirs, ids)

  known_ids <- unique(map$fastq_id)
  attribute <- function(df) {
    df$fastq_id <- vapply(df$filename, infer_fastq_id, character(1),
                          known_ids = known_ids)
    joined <- dplyr::left_join(df, map, by = "fastq_id")
    joined
  }
  stats <- attribute(fq$stats)
  modules <- attribute(fq$modules)

  unmatched <- unique(stats$filename[is.na(stats$run_id)])
  if (length(unmatched)) {
    warning(length(unmatched),
            " FastQC file(s) matched no fastq_id in any config.csv (e.g. '",
            unmatched[1], "').", call. = FALSE)
  }
  list(stats = stats, modules = modules)
}

#' Locate a feature-barcode HDF5 matrix for a run
#'
#' @param info Output of [detect_pipeline()].
#' @param which Either `"filtered"` (per-sample cell matrix) or `"raw"`
#'   (unfiltered GEM-well matrix).
#' @param sample_id Optional sample to select for `which = "filtered"`.
#' @return An absolute path, or `NA_character_` if none is found.
#' @keywords internal
#' @noRd
find_matrix_h5 <- function(info, which = c("filtered", "raw"), sample_id = NULL) {
  which <- match.arg(which)
  if (which == "raw") {
    candidates <- c(
      file.path(info$outs, "multi", "count", "raw_feature_bc_matrix.h5"),
      file.path(info$outs, "raw_feature_bc_matrix.h5")
    )
    hit <- candidates[file.exists(candidates)]
    return(if (length(hit)) normalizePath(hit[1]) else NA_character_)
  }

  samples <- info$samples
  if (!is.null(sample_id)) {
    samples <- samples[samples$sample_id == sample_id, , drop = FALSE]
    if (nrow(samples) == 0) {
      stop("Sample '", sample_id, "' not found in run.", call. = FALSE)
    }
  }
  for (i in seq_len(nrow(samples))) {
    candidates <- c(
      file.path(samples$count_dir[i], "sample_filtered_feature_bc_matrix.h5"),
      file.path(info$outs, "filtered_feature_bc_matrix.h5")
    )
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) return(normalizePath(hit[1]))
  }
  NA_character_
}

#' Read Cell Ranger version and chemistry from an h5 matrix
#' @keywords internal
#' @noRd
read_h5_run_attrs <- function(path) {
  rlang::check_installed("hdf5r")
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  attr_names <- hdf5r::h5attr_names(f)
  get_attr <- function(name) {
    if (name %in% attr_names) {
      as.character(hdf5r::h5attr(f, name))[1]
    } else {
      NA_character_
    }
  }
  list(
    software_version = get_attr("software_version"),
    chemistry_description = get_attr("chemistry_description")
  )
}

#' Read a 10x feature-barcode HDF5 matrix
#'
#' @param path Path to a `*_feature_bc_matrix.h5` file.
#' @return A list with `mat` (a features x barcodes `Matrix::dgCMatrix`),
#'   `features` (a tibble of `id`, `name`, `feature_type`, `genome`) and
#'   `barcodes`.
#' @keywords internal
#' @noRd
read_10x_h5 <- function(path) {
  rlang::check_installed(c("hdf5r", "Matrix"))
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  g <- f[["matrix"]]

  data <- g[["data"]]$read()
  indices <- g[["indices"]]$read()
  indptr <- g[["indptr"]]$read()
  shape <- g[["shape"]]$read()
  barcodes <- as.character(g[["barcodes"]]$read())

  ft <- g[["features"]]
  feature_type <- as.character(ft[["feature_type"]]$read())
  feature_name <- as.character(ft[["name"]]$read())
  feature_id <- as.character(ft[["id"]]$read())
  # `genome` is present for gene-expression references but may be absent.
  genome <- if ("genome" %in% names(ft)) {
    as.character(ft[["genome"]]$read())
  } else {
    rep(NA_character_, length(feature_id))
  }

  mat <- Matrix::sparseMatrix(
    i = as.numeric(indices) + 1,
    p = as.numeric(indptr),
    x = as.numeric(data),
    dims = as.integer(shape),
    dimnames = list(feature_id, barcodes)
  )

  list(
    mat = mat,
    features = tibble::tibble(
      id = feature_id,
      name = feature_name,
      feature_type = feature_type,
      genome = genome
    ),
    barcodes = barcodes
  )
}

#' Read a Cell Ranger feature-barcode matrix into a sparse matrix
#'
#' @param dir Path to a run directory (root or its `outs/`).
#' @param feature_type Optional feature type to keep (e.g.
#'   `"Gene Expression"`). If `NULL`, all features are returned.
#' @param sample_id For multi runs, which sample's filtered matrix to read.
#'   Defaults to the first sample.
#' @param which `"filtered"` (default) or `"raw"`.
#'
#' @return A features x barcodes `Matrix::dgCMatrix`.
#' @export
read_cr_matrix <- function(dir, feature_type = NULL, sample_id = NULL,
                           which = c("filtered", "raw")) {
  which <- match.arg(which)
  info <- detect_pipeline(dir)
  path <- find_matrix_h5(info, which = which, sample_id = sample_id)
  if (is.na(path)) {
    stop("No ", which, " feature-barcode .h5 found for run '", dir, "'.",
         call. = FALSE)
  }
  parsed <- read_10x_h5(path)
  if (!is.null(feature_type)) {
    keep <- parsed$features$feature_type %in% feature_type
    if (!any(keep)) {
      stop("No features of type '", paste(feature_type, collapse = ", "),
           "' in matrix.", call. = FALSE)
    }
    return(parsed$mat[keep, , drop = FALSE])
  }
  parsed$mat
}

#' Read a Cell Ranger matrix together with its feature metadata
#'
#' Like [read_cr_matrix()] but keeps the feature tibble (needed for
#' gene-symbol matching, e.g. mitochondrial/ribosomal fractions).
#'
#' @inheritParams read_cr_matrix
#' @return A list with `mat` (features x barcodes) and `features` (a tibble of
#'   `id`, `name`, `feature_type`, `genome`), both subset to `feature_type`.
#' @keywords internal
#' @noRd
read_cr_matrix_features <- function(dir, feature_type = NULL, sample_id = NULL,
                                    which = c("filtered", "raw")) {
  which <- match.arg(which)
  info <- detect_pipeline(dir)
  path <- find_matrix_h5(info, which = which, sample_id = sample_id)
  if (is.na(path)) {
    stop("No ", which, " feature-barcode .h5 found for run '", dir, "'.",
         call. = FALSE)
  }
  parsed <- read_10x_h5(path)
  if (!is.null(feature_type)) {
    keep <- parsed$features$feature_type %in% feature_type
    if (!any(keep)) {
      stop("No features of type '", paste(feature_type, collapse = ", "),
           "' in matrix.", call. = FALSE)
    }
    parsed$mat <- parsed$mat[keep, , drop = FALSE]
    parsed$features <- parsed$features[keep, , drop = FALSE]
  }
  list(mat = parsed$mat, features = parsed$features)
}

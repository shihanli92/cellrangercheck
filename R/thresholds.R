#' Default QC thresholds for Cell Ranger metrics
#'
#' A starting set of pass/warn/fail thresholds keyed by library type and metric.
#' `direction` is `"min"` when higher values are better (value below `warn`
#' warns, below `fail` fails) or `"max"` when lower is better. Percentage
#' metrics are on the 0-100 scale, matching [parse_cellranger()].
#'
#' Edit the returned tibble (add, drop or retune rows) and pass it to
#' [flag_qc()] to customise. `category`/`grouped_by` disambiguate which row of
#' the long metrics table a threshold applies to; `NA` in `grouped_by` matches
#' the per-sample `Cells` metrics.
#'
#' @return A tibble of thresholds with columns `library_type`, `category`,
#'   `grouped_by`, `metric`, `direction`, `warn`, `fail`.
#' @export
#' @examples
#' default_thresholds()
default_thresholds <- function() {
  tibble::tribble(
    ~library_type,      ~category, ~grouped_by,           ~metric,                               ~direction, ~warn,  ~fail,
    "Gene Expression",  "Cells",   NA,                    "Median genes per cell",               "min",      1000,   500,
    "Gene Expression",  "Cells",   NA,                    "Median UMI counts per cell",          "min",      2000,   1000,
    "Gene Expression",  "Cells",   NA,                    "Confidently mapped reads in cells",   "min",      70,     50,
    "Gene Expression",  "Library", "Physical library ID", "Valid barcodes",                      "min",      90,     75,
    "Gene Expression",  "Library", "Physical library ID", "Confidently mapped to transcriptome", "min",      60,     40,
    "Gene Expression",  "Library", "Physical library ID", "Q30 RNA read",                        "min",      80,     65,
    "Antibody Capture", "Library", "Physical library ID", "Valid barcodes",                      "min",      90,     75,
    "Antibody Capture", "Library", "Physical library ID", "Fraction antibody reads usable",      "min",      20,     10,
    "VDJ T",            "Cells",   NA,                    "Cells With Productive V-J Spanning Pair", "min",  30,     10,
    "VDJ T",            "Cells",   NA,                    "Fraction Reads in Cells",             "min",      50,     30,
    "Gene Expression",  "Derived", NA,                    "Median % mitochondrial",              "max",      10,     20,
    "Gene Expression",  "Derived", NA,                    "% cells passing filter",              "min",      70,     50,
    "Gene Expression",  "Derived", NA,                    "% ambient RNA",                       "max",      20,     40
  )
}

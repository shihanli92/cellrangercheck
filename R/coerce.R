#' Coerce a Cell Ranger metric value string to a number and unit
#'
#' Cell Ranger writes metric values as strings such as `"86.3%"`, `"47,538"`,
#' `"2.0"` or `"412.07"`. This strips thousands separators and a trailing
#' percent sign and classifies the unit.
#'
#' @param x Character vector of raw metric values.
#'
#' @return A list with two equal-length elements: `value` (numeric) and
#'   `unit` (character, one of `"percent"`, `"count"` or `"number"`). Percent
#'   values are returned on their original 0-100 scale (e.g. `86.3`), not as a
#'   fraction; the `unit` marks them so consumers can rescale if desired.
#' @keywords internal
#' @noRd
coerce_metric_value <- function(x) {
  x <- as.character(x)
  is_pct <- grepl("%", x, fixed = TRUE)
  # A value is a "count" if it carried a thousands separator (integers like
  # cell counts and read counts); everything else is a plain number.
  had_comma <- grepl(",", x, fixed = TRUE)
  cleaned <- gsub("[,%]", "", x)
  cleaned <- trimws(cleaned)
  cleaned[cleaned == ""] <- NA_character_
  value <- suppressWarnings(as.numeric(cleaned))
  unit <- ifelse(is_pct, "percent", ifelse(had_comma, "count", "number"))
  unit[is.na(value)] <- NA_character_
  list(value = value, unit = unit)
}

#' Turn a metric label into a snake_case token
#'
#' @param x Character vector of human-readable metric names.
#' @return Character vector of snake_case tokens.
#' @keywords internal
#' @noRd
snakecase <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

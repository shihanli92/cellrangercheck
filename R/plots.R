#' Tidy per-assay metric table for plotting
#'
#' Reshapes the long metrics table into one row per (library, metric) for use in
#' the metric-explorer plots: per-sample cell metrics and per-library sequencing
#' metrics, labelled `[cell]` / `[seq]` so same-named metrics stay distinct.
#'
#' @param long A tidy tibble from [parse_cellranger_runs()].
#'
#' @return A tibble with `run_id`, `sample_id`, `library_type`, `library`,
#'   `metric_label`, `value`, `unit`.
#' @export
assay_metric_data <- function(long) {
  keep <-
    long$category %in% c("Cells", "Derived") |
    (long$category %in% "Library" & long$grouped_by %in% "Physical library ID") |
    is.na(long$category)
  d <- long[keep & !is.na(long$value), , drop = FALSE]

  tag <- ifelse(is.na(d$category), "",
                ifelse(d$category == "Cells", " [cell]",
                       ifelse(d$category == "Derived", " [derived]", " [seq]")))
  d$metric_label <- paste0(d$metric, tag)
  d$library <- ifelse(
    is.na(d$sample_id) | d$run_id == d$sample_id,
    d$run_id,
    paste(d$run_id, d$sample_id, sep = "/")
  )
  tibble::as_tibble(
    d[c("run_id", "sample_id", "library_type", "library",
        "metric_label", "value", "unit")]
  )
}

#' Interactive library-vs-metric bar plot for one assay
#'
#' Builds a plotly bar chart with a dropdown to choose which numeric metric to
#' show: x = library, y = the selected metric. Requires the suggested package
#' `plotly`.
#'
#' @param mdat Output of [assay_metric_data()].
#' @param library_type Which assay (library type) to plot.
#'
#' @return A plotly htmlwidget.
#' @export
assay_metric_plot <- function(mdat, library_type) {
  rlang::check_installed("plotly")
  d <- mdat[mdat$library_type == library_type, , drop = FALSE]
  if (nrow(d) == 0) {
    return(plotly::plotly_empty(type = "scatter", mode = "markers"))
  }
  metrics <- sort(unique(d$metric_label))
  libs <- sort(unique(d$library))

  p <- plotly::plot_ly()
  for (i in seq_along(metrics)) {
    dm <- d[d$metric_label == metrics[i], , drop = FALSE]
    dm <- dm[match(libs, dm$library), , drop = FALSE]
    p <- plotly::add_bars(
      p, x = libs, y = dm$value, name = metrics[i],
      visible = (i == 1),
      hovertemplate = "%{x}<br>%{y}<extra></extra>"
    )
  }
  buttons <- lapply(seq_along(metrics), function(i) {
    vis <- as.list(rep(FALSE, length(metrics)))
    vis[[i]] <- TRUE
    list(method = "update", label = metrics[i],
         args = list(list(visible = vis),
                     list(yaxis = list(title = metrics[i]))))
  })
  plotly::layout(
    p,
    title = library_type,
    xaxis = list(title = "library"),
    yaxis = list(title = metrics[1]),
    showlegend = FALSE,
    updatemenus = list(list(
      active = 0, buttons = buttons, direction = "down",
      x = 1.02, xanchor = "left", y = 1, yanchor = "top",
      pad = list(r = 5)
    )),
    margin = list(r = 160)
  )
}

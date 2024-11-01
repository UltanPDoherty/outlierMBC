#' @title Find the gross outliers.
#'
#' @inheritParams ombc_gmm
#' @param multiplier Factor by which the `max_out`th kNN Distance is multiplied
#'                   to get the gross outlier threshold.
#' @param k_neighbours Number of neighbours for dbscan::kNNdist.
#' @param manual_gross_threshold Optional preset number of gross outliers.
#' @param scale Logical
#'
#' @return List:
#' * $choice: a numeric value indicating the elbow's location.
#' * $bool: a logical vector identifying the gross outliers.
#' * $plot
#'
#' @export
find_gross <- function(
    x,
    max_out,
    multiplier = 3,
    k_neighbours = floor(nrow(x) / 100),
    manual_gross_threshold = NULL,
    scale = TRUE) {
  outlier_number <- seq_len(2 * max_out)

  if (scale) {
    x <- scale(x)
  }

  x_knndist <- dbscan::kNNdist(x, k_neighbours)
  knndist_sort <- -sort(-x_knndist)[outlier_number]

  knndist_maxout <- knndist_sort[max_out]

  if (is.null(manual_gross_threshold)) {
    gross_threshold <- multiplier * knndist_maxout
  } else {
    gross_threshold <- manual_gross_threshold
    multiplier <- NA
  }

  gross_bool <- x_knndist > gross_threshold
  gross_choice <- sum(gross_bool)

  gross <- NULL
  curve <- data.frame(
    outlier_number, knndist_sort,
    gross = outlier_number <= gross_choice
  ) |>
    ggplot2::ggplot(
      ggplot2::aes(x = outlier_number, y = knndist_sort, colour = gross)
    ) +
    ggplot2::geom_point(
      size = min(1, max(0.1, 100 / max_out)), show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = gross_threshold, colour = "#E69F00") +
    ggplot2::labs(
      x = "kNN Distance Order",
      y = paste0("kNN Distance (k = ", k_neighbours, ")"),
      colour = "Outlier Number Choices:"
    ) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_colour_manual(values = c("#000000", "#E69F00"))

  if (is.null(manual_gross_threshold)) {
    curve <- curve +
      ggplot2::geom_linerange(
        ymax = knndist_maxout, ymin = 0, x = max_out,
        linetype = "dashed", colour = "black", show.legend = FALSE
      ) +
      ggplot2::geom_linerange(
        y = knndist_maxout, xmin = 0, xmax = max_out,
        linetype = "dashed", colour = "black", show.legend = FALSE
      )
  }

  x_seq <- NULL
  scatter <- data.frame(x_seq = seq_len(nrow(x)), x_knndist, gross_bool) |>
    ggplot2::ggplot(ggplot2::aes(
      x = x_seq, y = x_knndist, colour = gross_bool
    )) +
    ggplot2::geom_point(
      size = min(1, max(0.1, 100 / max_out)), show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 3 * knndist_maxout, colour = "#E69F00") +
    ggplot2::labs(
      x = "Index",
      y = paste0("kNN Distance (k = ", k_neighbours, ")")
    ) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_colour_manual(values = c("#000000", "#E69F00"))

  if (is.null(manual_gross_threshold)) {
    scatter <- scatter +
      ggplot2::geom_hline(yintercept = knndist_maxout, linetype = "dashed")
  }

  plot <- ggpubr::annotate_figure(
    ggpubr::ggarrange(curve, scatter, nrow = 1, ncol = 2),
    top = paste0(
      "No. of Gross Outliers = ", gross_choice,
      " (max_out = ", max_out, ", multiplier = ", multiplier, ")"
    )
  ) + ggpubr::bgcolor("white")

  output <- list(
    gross_choice = gross_choice,
    gross_bool = gross_bool,
    plot = plot
  )

  return(output)
}

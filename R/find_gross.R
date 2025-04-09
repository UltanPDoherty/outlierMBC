#' @title Find gross outliers.
#'
#' @description
#' The distance of each observation to its \eqn{k^{th}}{k^th} nearest neighbour
#' is computed. We assume that the largest `max_out` kNN distances correspond to
#' potential outliers. We select the next largest kNN distance, outside of the
#' top `max_out`, as a benchmark value. We multiply this benchmark kNN distance
#' by `multiplier` to get the minimum threshold for our gross outliers. In other
#' words, a gross outlier must have a kNN distance at least `multiplier` times
#' greater than all of the observations which we do not consider to be potential
#' outliers.
#'
#' @inheritParams ombc_gmm
#' @param multiplier Multiplicative factor used to get gross outlier threshold.
#' @param k_neighbours Number of neighbours for dbscan::kNNdist.
#' @param manual_threshold Optional preset threshold.
#' @param scale Logical value controlling whether we apply `scale` to `x`.
#'
#' @return
#' `find_gross` returns a list with the following elements:
#' \describe{
#'   \item{`gross_choice`}{A numeric value indicating the elbow's location.}
#'   \item{`gross_bool`}{A logical vector identifying the gross outliers.}
#'   \item{`gross_curve`}{ggplot of the highest `2 * max_out` kNN distances in
#'                 decreasing order.}
#'   \item{`gross_scatter`}{ggplot of all kNN distances in index order.}
#' }
#'
#' @export
find_gross <- function(
    x,
    max_out,
    multiplier = 3,
    k_neighbours = floor(nrow(x) / 100),
    manual_threshold = NULL,
    scale = TRUE) {
  outlier_number <- seq_len(2 * max_out)

  if (scale) {
    x <- scale(x)
  }

  x_knndist <- dbscan::kNNdist(x, k_neighbours)
  knndist_sort <- -sort(-x_knndist)[outlier_number]

  knndist_benchmark <- knndist_sort[max_out + 1]

  if (is.null(manual_threshold)) {
    gross_threshold <- multiplier * knndist_benchmark
  } else {
    gross_threshold <- manual_threshold
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
      title = paste0(
        "No. of Gross Outliers = ", gross_choice,
        " (max_out = ", max_out, ", multiplier = ", multiplier, ")"
      ),
      subtitle = paste0(
        2 * max_out, " largest kNN distances (k = ", k_neighbours, ")"
      ),
      x = "kNN Distance Order",
      y = paste0("kNN Distance (k = ", k_neighbours, ")"),
      colour = "Outlier Number Choices:"
    ) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_colour_manual(values = c("#000000", "#E69F00"))

  if (is.null(manual_threshold)) {
    curve <- curve +
      ggplot2::geom_linerange(
        ymax = knndist_benchmark, ymin = 0, x = max_out + 1,
        linetype = "dashed", colour = "black", show.legend = FALSE
      ) +
      ggplot2::geom_linerange(
        y = knndist_benchmark, xmin = 0, xmax = max_out + 1,
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
    ggplot2::geom_hline(yintercept = gross_threshold, colour = "#E69F00") +
    ggplot2::labs(
      title = paste0(
        "No. of Gross Outliers = ", gross_choice,
        " (max_out = ", max_out, ", multiplier = ", multiplier, ")"
      ),
      subtitle = paste0(
        "All kNN distances (k = ", k_neighbours, ")"
      ),
      x = "Index",
      y = paste0("kNN Distance (k = ", k_neighbours, ")")
    ) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_colour_manual(values = c("#000000", "#E69F00"))

  if (is.null(manual_threshold)) {
    scatter <- scatter +
      ggplot2::geom_hline(yintercept = knndist_benchmark, linetype = "dashed")
  }

  list(
    gross_choice = gross_choice,
    gross_bool = gross_bool,
    gross_curve = curve,
    gross_scatter = scatter
  )
}

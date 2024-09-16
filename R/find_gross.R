#' @title Find the gross outliers.
#'
#' @inheritParams find_elbow
#' @inheritParams outcast_gmm
#' @param k_neighbours Number of neighbours for dbscan::glosh
#' @param choice Optional preset number of gross outliers.
#'
#' @return List:
#' * $choice: a numeric value indicating the elbow's location.
#' * $bool: a logical vector identifying the gross outliers.
#' * $plot
#'
#' @export
find_gross <- function(
    x, max_out, search_centre, k_neighbours = 10, choice = NULL) {
  x_glosh <- dbscan::glosh(x, k_neighbours)
  glosh_sort <- -sort(-x_glosh)[seq_len(max_out)]

  elbow <- find_elbow(glosh_sort, search_centre, FALSE)

  if (is.null(choice)) {
    choice <- elbow$choice
  }

  bool <- rank(-x_glosh) <= choice

  return(list(choice = choice, bool = bool, plot = elbow$plot))
}

#' Plot sorted dbscan::glosh values to determine the `search_centre` for
#' `find_gross`.
#'
#' @inheritParams find_gross
#'
#' @return A ggplot2 object
#'
#' @export
plot_gross <- function(x, max_out, k_neighbours = 10) {
  x_glosh <- dbscan::glosh(x, k_neighbours)

  outlier_seq <- seq_len(max_out)
  glosh_sort <- -sort(-x_glosh)[seq_len(max_out)]
  gg <- data.frame(outlier_seq, glosh_sort) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = glosh_sort)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Outlier Number",
      y = paste0("Sorted dbscan::glosh Values (k = ", k_neighbours, ")")
    )

  return(gg)
}

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
    x, search_centre,
    k_neighbours = floor(nrow(x) / 100),
    underestimate = 0.5,
    choice = NULL) {
  
  stopifnot(!is.null(search_centre))
  
  x_knndist <- dbscan::kNNdist(x, k_neighbours)
  knndist_sort <- -sort(-x_knndist)
  
  elbow <- find_elbow(knndist_sort, search_centre, TRUE)

  if (is.null(choice)) {
    choice <- elbow$choice
  }
  
  gross_choice <- floor(choice * underestimate)

  bool <- rank(-x_knndist) <= gross_choice
  
  plot <- elbow$plot + 
    ggplot2::geom_vline(xintercept = gross_choice, colour = "red") +
    ggplot2::labs(
      subtitle = paste0(
        "Gross choice = ", gross_choice
      ),
      y = paste0("kNN Distance (k = ", k_neighbours, ")")
    )

  return(list(choice = gross_choice, bool = bool, plot = plot))
}

#' Plot sorted dbscan::glosh values to determine the `search_centre` for
#' `find_gross`.
#'
#' @inheritParams find_gross
#'
#' @return A ggplot2 object
#'
#' @export
plot_gross <- function(x, k_neighbours = floor(nrow(x) / 100)) {

  x_knndist <- dbscan::kNNdist(x, k_neighbours)
  
  outlier_seq <- seq_len(nrow(x))
  knndist_sort <- -sort(-x_knndist)
  gg <- data.frame(outlier_seq, knndist_sort) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = knndist_sort)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Outlier Number",
      y = paste0("kNN Distance (k = ", k_neighbours, ")")
    )

  return(gg)
}

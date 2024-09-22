#' @title Find the gross outliers.
#'
#' @inheritParams ombc1_gmm
#' @param k_neighbours Number of neighbours for dbscan::kNNdist.
#' @param underestimate Factor by which to multiply the elbow estimate of the
#'                      number of outliers to get the number of gross outliers.
#' @param search_centre Optional centre of elbow search interval.
#' @param choice Optional preset number of gross outliers.
#'
#' @return List:
#' * $choice: a numeric value indicating the elbow's location.
#' * $bool: a logical vector identifying the gross outliers.
#' * $plot
#'
#' @export
find_gross <- function(
    x, max_out,
    k_neighbours = floor(nrow(x) / 100),
    underestimate = 0.5,
    search_centre = NULL,
    choice = NULL) {

  outlier_number <- seq_len(2 * max_out)

  x_knndist <- dbscan::kNNdist(x, k_neighbours)
  knndist_sort <- -sort(-x_knndist)[outlier_number]

  if (is.null(search_centre)) {
     cpts <- changepoint::cpt.meanvar(diff(knndist_sort))@cpts
     stopifnot(
       "No suitable search_centre found.\n" = length(cpts) > 1
      )
  }

  elbow <- find_elbow(knndist_sort, search_centre, TRUE)
  elbow_choice <- elbow$choice

  gross_choice <- floor(elbow_choice * underestimate)

  if (!is.null(choice)) {
    gross_choice <- choice
  }

  bool <- rank(-x_knndist) <= gross_choice

  search_interval <- elbow$search_interval
  cpop_fitted <- cpop::fitted(elbow$cpop_out)
  gg <- data.frame(outlier_number, knndist_sort) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_number, y = knndist_sort)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = search_interval, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = elbow_choice) +
    ggplot2::geom_vline(xintercept = gross_choice, colour = "red") +
    ggplot2::geom_abline(
      slope = cpop_fitted$gradient[1], intercept = cpop_fitted$intercept[1],
      linetype = "dotted"
    ) +
    ggplot2::geom_abline(
      slope = cpop_fitted$gradient[2], intercept = cpop_fitted$intercept[2],
      linetype = "dotted"
    ) +
    ggplot2::labs(
      title = paste0(
        "Chosen number of gross outliers = ", gross_choice
      ),
      subtitle = paste0(
        "Elbow choice = ", elbow_choice,
        " (search interval = [", search_interval[1],
        ", ", search_interval[2], "])"
      ),
      x = "Outlier Number",
      y = paste0("kNN Distance (k = ", k_neighbours, ")")
    )

  return(list(choice = gross_choice, bool = bool, plot = gg))
}

#' Plot kNN distance values to determine the `search_centre` for
#' `find_gross`.
#'
#' @inheritParams find_gross
#'
#' @return A ggplot2 object
#'
#' @export
plot_gross <- function(x, max_out, k_neighbours = floor(nrow(x) / 100)) {
  x_knndist <- dbscan::kNNdist(x, k_neighbours)

  outlier_seq <- seq_len(2 * max_out)
  knndist_sort <- -sort(-x_knndist)[outlier_seq]
  gg <- data.frame(outlier_seq, knndist_sort) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = knndist_sort)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Outlier Number",
      y = paste0("kNN Distance (k = ", k_neighbours, ")")
    )

  return(gg)
}

find_elbow <- function(y, search_centre = NULL, concave = TRUE) {
  y_len <- length(y)
  y_seq <- seq_len(y_len)

  if (is.null(search_centre)) {
    linmod <- stats::lm(y ~ y_seq)
    if (concave) {
      search_centre <- which.min(linmod$residuals)
    } else {
      search_centre <- which.max(linmod$residuals)
    }
  }

  search_radius <- min(c(
    floor((y_len - search_centre - 1) / 3),
    floor((search_centre - 2) / 3)
  ))

  if (search_radius == 0 && search_centre < 5) {
    warning("search_radius == 0. Try search_centre > 4.\n")
  } else if (search_radius == 0 && search_centre > (y_len - 4)) {
    warning(paste0(
      "search_radius == 0. Try search_centre < ", y_len - 3, ".\n"
    ))
  }

  upper <- search_centre + search_radius
  lower <- search_centre - search_radius
  search_interval <- c(lower, upper)

  width <- floor(y_len / 10)
  left <- right <- integer(y_len)
  for (i in y_seq) {
    left[i] <- min(max(c(i - width, 1)), y_len - 2 * width)
    right[i] <- max(min(c(i + width, y_len)), 1 + 2 * width)
  }

  cpop_out <- suppressMessages(cpop::cpop(
    y, y_seq,
    grid = lower:upper, minseglen = 2 * search_radius + 1
  ))

  stopifnot(length(cpop_out@changepoints) == 3)

  choice <- cpop_out@changepoints[2]

  cat("search centre = ", search_centre)
  cat(", search radius = ", search_radius)
  cat(", search interval = ", search_interval)
  cat(", search choice = ", choice, "\n")

  return(list(
    choice = choice,
    search_interval = search_interval,
    cpop_out = cpop_out
  ))
}

# ==============================================================================

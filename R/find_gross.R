#' @title Find the gross outliers.
#'
#' @inheritParams ombc_gmm
#' @param k_neighbours Number of neighbours for dbscan::kNNdist.
#' @param gross_prop Factor by which to multiply the elbow estimate of the
#'                   number of outliers to get the number of gross outliers.
#' @param search_centre Optional centre of elbow search interval.
#' @param manual_choice Optional preset number of gross outliers.
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
    gross_prop = max(0.6, min(0.8, 0.9 - (30 / elbow_choice))),
    search_centre = NULL,
    manual_choice = NULL) {
  outlier_number <- seq_len(2 * max_out)

  x_knndist <- dbscan::kNNdist(x, k_neighbours)
  knndist_sort <- -sort(-x_knndist)[outlier_number]

  if (is.null(search_centre)) {
    cpts <- changepoint::cpt.meanvar(diff(knndist_sort))@cpts
    stopifnot("No suitable search_centre found.\n" = length(cpts) > 1)
    search_centre <- cpts[1]
  }

  elbow <- find_elbow(knndist_sort, search_centre, TRUE)
  elbow_choice <- elbow$choice

  gross_choice <- floor(elbow_choice * gross_prop)

  if (!is.null(manual_choice)) {
    gross_choice <- manual_choice

    cat(paste0("manual choice = ", manual_choice, ".\n"))
  } else {
    cat(paste0(
      "elbow = ", elbow_choice,
      ", gross_prop = ", round(gross_prop, 2),
      ", no. of gross outliers selected = ", gross_choice, ".\n"
    ))
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
        ", ", search_interval[2],
        "]), gross proportion = ", round(gross_prop, 2)
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

  cat(paste0(
    "elbow search: centre = ", search_centre,
    ", radius = ", search_radius,
    ", interval = (", search_interval[1], ", ", search_interval[2], ")",
    ", choice = ", choice, ".\n"
  ))

  return(list(
    choice = choice,
    search_interval = search_interval,
    cpop_out = cpop_out
  ))
}

# ==============================================================================

#' @title Find the gross outliers.
#'
#' @inheritParams ombc_gmm
#' @param k_neighbours Number of neighbours for dbscan::kNNdist.
#' @param gross_prop Factor by which to multiply the elbow estimate of the
#'                   number of outliers to get the number of gross outliers.
#' @param search_centre Optional centre of elbow search interval.
#' @param manual_choice Optional preset number of gross outliers.
#'
#' @return List:
#' * $choice: a numeric value indicating the elbow's location.
#' * $bool: a logical vector identifying the gross outliers.
#' * $plot
#'
#' @export
find_gross2 <- function(
    x, max_out,
    k_neighbours = floor(nrow(x) / 100),
    gross_prop = max(0.6, min(0.8, 0.9 - (30 / elbow_choice))),
    candidate,
    manual_choice = NULL) {
  outlier_number <- seq_len(2 * max_out)

  x_knndist <- dbscan::kNNdist(x, k_neighbours)
  knndist_sort <- -sort(-x_knndist)[outlier_number]

  candidates <- c()

  cpts_meanvar <- changepoint::cpt.meanvar(knndist_sort)@cpts
  if (length(cpts_meanvar) == 2) {
    candidates <- append(candidates, cpts_meanvar[1])
    candidates <- append(candidates, cpts_meanvar[1] + 1)
  }

  cpts_diff_meanvar <- changepoint::cpt.meanvar(diff(knndist_sort))@cpts
  if (length(cpts_diff_meanvar) == 2) {
    candidates <- append(candidates, cpts_diff_meanvar[1])
    candidates <- append(candidates, cpts_diff_meanvar[1] + 1)
  }

  stopifnot(
    "No suitable candidates found.\n" = length(candidates) > 0
  )

  cat(paste0(
    "Initial elbow candidates: ", paste0(candidates, collapse = ", "), ".\n"
  ))

  candidates <- unique(append(
    candidates,
    round(seq(
      max(1, min(candidates) - 10),
      min(max_out, max(candidates) + 10),
      length.out = 20
    ))
  ))

  cat(paste0(
    "All elbow candidates:\n\t",
    paste0(
      sort(candidates)[seq(1, 1 + (length(candidates) %/% 2))],
      collapse = ", "
    ),
    ",\n\t",
    paste0(
      sort(candidates)[seq(2 + (length(candidates) %/% 2), length(candidates))],
      collapse = ", "
    ),
    ".\n"
  ))

  tests <- lapply(
    candidates, function(t) lm_test(knndist_sort, outlier_number, t)
  )

  test_scores <- vapply(tests, function(x) x$rss, double(1L))
  which_best <- which.min(test_scores)

  best_test <- tests[[which_best]]
  elbow_choice <- candidates[which_best]

  gross_choice <- floor(elbow_choice * gross_prop)

  if (!is.null(manual_choice)) {
    gross_choice <- manual_choice

    cat(paste0("manual choice = ", manual_choice, ".\n"))
  } else {
    cat(paste0(
      "Elbow choice = ", elbow_choice,
      ", gross proportion = ", round(gross_prop, 2),
      ", no. of gross outliers selected = ", gross_choice, ".\n"
    ))
  }

  elbow_bool <- rank(-x_knndist) <= elbow_choice
  gross_bool <- rank(-x_knndist) <= gross_choice

  choice_df <- data.frame("elbow" = elbow_choice, "gross" = gross_choice)
  gg <- data.frame(outlier_number, knndist_sort) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_number, y = knndist_sort)) +
    ggplot2::geom_point(size = min(1, max(0.1, 100 / max_out))) +
    ggplot2::geom_vline(
      data = choice_df,
      ggplot2::aes(xintercept = elbow, colour = "elbow"), linetype = "dashed",
      linewidth = 0.75
    ) +
    ggplot2::geom_vline(
      data = choice_df,
      ggplot2::aes(xintercept = gross, colour = "gross"), linetype = "dashed",
      linewidth = 0.75
    ) +
    ggplot2::geom_segment(
      x = 0, y = best_test$coeff1[1],
      xend = elbow_choice,
      yend = best_test$coeff1[1] + elbow_choice * best_test$coeff1[2],
      linewidth = 0.75, colour = "#F0E442"
    ) +
    ggplot2::geom_segment(
      x = elbow_choice + 1,
      y = best_test$coeff2[1] + elbow_choice * best_test$coeff2[2],
      xend = 2 * max_out,
      yend = best_test$coeff2[1] + 2 * max_out * best_test$coeff2[2],
      linewidth = 0.75, colour = "#F0E442"
    ) +
    ggplot2::labs(
      title = paste0(
        "Chosen number of gross outliers = ", gross_choice
      ),
      subtitle = paste0(
        "Elbow choice = ", elbow_choice,
        ", gross proportion = ", round(gross_prop, 2)
      ),
      x = "Outlier Number",
      y = paste0("kNN Distance (k = ", k_neighbours, ")"),
      colour = "Outlier Number Choices:"
    ) +
    ggplot2::scale_colour_manual(
      values = c(elbow = "#56B4E9", gross = "#E69F00")
    ) +
    ggplot2::theme(legend.position = "bottom")

  return(list(
    gross_choice = gross_choice,
    gross_bool = gross_bool,
    elbow_choice = elbow_choice,
    elbow_bool = elbow_bool,
    plot = gg
  ))
}

lm_test <- function(y, x, split) {
  y1 <- y[seq(1, split)]
  x1 <- x[seq(1, split)]

  y2 <- y[seq(split+1, length(y))]
  x2 <- x[seq(split+1, length(x))]

  lm1 <- lm(y1 ~ x1)
  lm2 <- lm(y2 ~ x2)

  rss1 <- sum(lm1$residuals^2)
  rss2 <- sum(lm2$residuals^2)

  return(list(
    rss = rss1 + rss2,
    coeff1 = lm1$coefficients,
    coeff2 = lm2$coefficients
  ))
}

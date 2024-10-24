#' @title Find the gross outliers.
#'
#' @inheritParams ombc_gmm
#' @param k_neighbours Number of neighbours for dbscan::kNNdist.
#' @param manual_gross_choice Optional preset number of gross outliers.
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
    manual_gross_choice = NULL) {
  outlier_number <- seq_len(2 * max_out)

  x_knndist <- dbscan::kNNdist(x, k_neighbours)
  knndist_sort <- -sort(-x_knndist)[outlier_number]

  candidates <- round(seq(2, max_out - 2, by = 1))

  tests <- lapply(
    candidates, function(t) lm_test(knndist_sort, outlier_number, t)
  )

  test_scores <- vapply(tests, function(x) x$rss, double(1L))

  if (all(is.infinite(test_scores))) {
    cat("No gross outliers identified.\n")

    elbow_choice <- 0
  } else {
    which_best <- which.min(test_scores)

    best_test <- tests[[which_best]]
    elbow_choice <- candidates[which_best]
  }

  gross_prop <- max(0.6, min(0.9, 0.95 - (35 / elbow_choice)))
  gross_choice <- floor(elbow_choice * gross_prop)

  if (!is.null(manual_gross_choice)) {
    gross_choice <- manual_gross_choice

    cat(paste0("manual gross choice = ", manual_gross_choice, ".\n"))
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
  elbow <- gross <- NULL
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

  if (!all(is.infinite(test_scores))) {
    gg <- gg +
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
      )
  }

  output <- list(
    gross_choice = gross_choice,
    gross_bool = gross_bool,
    elbow_choice = elbow_choice,
    elbow_bool = elbow_bool,
    plot = gg
  )

  return(output)
}

lm_test <- function(y, x, split) {
  df1 <- data.frame("x1" = x[seq(1, split)], "y1" = y[seq(1, split)])
  df2 <- data.frame("x2" = x[-seq(1, split)], "y2" = y[-seq(1, split)])

  lm0 <- lm(y ~ x)
  rss0 <- sum(lm0$residuals^2)

  lm1 <- lm(y1 ~ x1, data = df1)
  rss1 <- sum(lm1$residuals^2)

  lm2 <- lm(y2 ~ x2, data = df2)
  rss2 <- sum(lm2$residuals^2)

  slope_ratio <- lm1$coefficients[2] / lm2$coefficients[2]
  if (slope_ratio > 10) {
    rss <- rss1 + rss2
  } else {
    rss <- Inf
  }

  rss_ratio <- rss0 / rss

  if (rss_ratio < 10) {
    rss <- Inf
  }

  return(list(
    rss = rss,
    coeff1 = lm1$coefficients,
    coeff2 = lm2$coefficients
  ))
}

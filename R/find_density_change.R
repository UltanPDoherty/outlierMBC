#' @export
find_density_change <- function(y, search_centre) {
  y_len <- length(y)
  outliers_removed <- seq_len(y_len)
  
  search_radius <- min(c(
    floor((y_len - search_centre - 1) / 3),
    floor((search_centre - 2) / 3)
  ))
  
  if (search_radius == 0 & search_centre < 5) {
    warning("search_radius == 0. Try search_centre > 4.\n")
  } else if (search_radius == 0 & search_centre > (y_len - 4)) {
    warning(paste0("search_radius == 0. Try search_centre < ", y_len - 3, ".\n"))
  }
  
  upper <- search_centre + search_radius
  lower <- search_centre - search_radius
  search_interval <- c(lower, upper)
  
  width <- floor(y_len / 10)
  left <- right <- integer(y_len)
  sd_vec <- double(y_len)
  for (i in outliers_removed) {
    left[i] <- min(max(c(i - width, 1)), y_len - 2 * width)
    right[i] <- max(min(c(i + width, y_len)), 1 + 2 * width)
    sd_vec[i] <- sqrt(mean(diff(diff(y[left[i]:right[i]]))^2) / 6)
  }
  
  cpop_out <- suppressMessages(cpop::cpop(
    y, outliers_removed,
    grid = lower:upper, minseglen = 2 * search_radius + 1,
    sd = sd_vec
  ))
  
  stopifnot(length(cpop_out@changepoints) == 3)
  
  choice <- cpop_out@changepoints[2]
  
  cat("search centre = ", search_centre)
  cat(", search radius = ", search_radius)
  cat(", search interval = ", search_interval)
  cat(", search choice = ", cpop_out@changepoints[2], "\n")
  
  cpop_fitted <- cpop::fitted(cpop_out)
  
  gg <- cbind(x = outliers_removed, y = y) |>
    ggplot2::ggplot(ggplot2::aes(x = outliers_removed, y = y)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = search_interval, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = choice) +
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
        "choice = ", choice,
        " (search_interval = [", search_interval[1],
        ", ", search_interval[2], "])"
      )
    )
  
  return(list(
    choice = choice,
    search_interval = search_interval,
    cpop_out = cpop_out,
    plot = gg
  ))
}

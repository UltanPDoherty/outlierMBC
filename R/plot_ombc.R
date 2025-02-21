#' Plot the outlier number selection curve.
#'
#' @param ombc_out Output from ombc_gmm.
#'
#' @returns A gg object.
#' @export
plot_curve <- function(ombc_out) {
  gross_num <- sum(ombc_out$gross_outs)
  max_out <- max(ombc_out$outlier_rank) - 1
  outlier_num <- ombc_out$outlier_num
  distrib_diff_vec <- ombc_out$distrib_diff_vec

  outlier_seq <- seq(gross_num, max_out)
  point_size <- 1 - min(0.9, max(0, -0.1 + max_out / 250))

  dd <- minimum <- NULL
  curve_df <- data.frame(
    "outlier_seq" = outlier_seq,
    "minimum" = as.integer(outlier_num),
    "dd" = distrib_diff_vec
  )
  curve <- curve_df |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = dd)) +
    ggplot2::geom_line(
      ggplot2::aes(colour = "dd"),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = "dd"),
      size = point_size, show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = minimum, colour = "minimum"),
      linetype = "solid", linewidth = 0.75, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(dd = "#000000", minimum = "#CC79A7")
    ) +
    ggplot2::labs(
      title = paste0("outlierMBC: Number of Outliers = ", outlier_num),
      x = "Outlier Number",
      y = "Dissimilarity",
      colour = ""
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(outlier_seq)) +
    ggplot2::expand_limits(y = 0)

  return(curve)
}

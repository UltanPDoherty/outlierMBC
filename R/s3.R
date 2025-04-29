#' Constructor for `"outliermbc_gmm"` S3 class.
#'
#' @param x List.
#'
#' @returns "outliermbc_gmm" S3 object.
new_outliermbc_gmm <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class = c("outliermbc_gmm", "list"))
}

# ------------------------------------------------------------------------------

#' Validator for `"outliermbc_gmm"` S3 class.
#'
#' @param x List.
#'
#' @returns NULL
validate_outliermbc_gmm <- function(x) {
  values <- unclass(x)

  obs_num <- lapply(
    values[c("labels", "outlier_bool", "outlier_rank", "gross_outs")],
    length
  )
  if (length(unique(obs_num)) != 1) {
    stop(
      paste0(
        "labels, outlier_bool, outlier_rank, and gross_outs should have the ",
        "same length."
      ),
      call. = FALSE
    )
  }

  iter_num <- c(lapply(
    values[c("loglike", "removal_dens", "distrib_diff_vec")],
    length
  ), nrow(values$distrib_diff_mat))
  if (length(unique(iter_num)) != 1) {
    stop(
      paste0(
        "loglike, removal_dens, and distrib_diff_vec should have the same ",
        "length. This should also equal the number of rows in distrib_diff_mat."
      ),
      call. = FALSE
    )
  }
}

# ------------------------------------------------------------------------------

#' plot method for `"outliermbc_gmm"` S3 class.
#'
#' @param x List
#' @param backtrack Logical
#' @param ... Other arguments
#'
#' @returns A ggplot
#' @export
plot.outliermbc_gmm <- function(x, backtrack = FALSE, ...) {
  if (!backtrack) {
    plot_curve(x)
  } else {
    plot_backtrack(x)
  }
}

# ------------------------------------------------------------------------------

#' print method for `"outliermbc_gmm"` S3 class.
#'
#' @param x List
#' @param backtrack Logical
#' @param ... Other arguments
#' @inheritParams backtrack
#'
#' @returns A ggplot
#' @export
print.outliermbc_gmm <- function(
    x, backtrack = FALSE, max_total_rise = 0.1, max_step_rise = 0.05, ...) {
  obs_num <- length(x$labels)
  max_out <- x$call$max_out
  gross_num <- sum(x$gross_outs)
  outlier_num <- x$outlier_num

  cat("Starting number of data points:\t", obs_num, "\n")
  cat("Maximum number of outliers:\t", max_out, "\n")
  cat("Number of gross outliers:\t", gross_num, "\n")

  if (!backtrack) {
    cat(
      "Final number of outliers:\t", outlier_num, "(minimum dissimilarity)\n"
    )
  } else {
    backtrack_out <- backtrack(
      x$distrib_diff_vec, max_total_rise, max_step_rise
    )
    backtrack_num <- backtrack_out$backtrack$ind - 1 + gross_num

    cat(
      "Provisional number of outliers:\t", outlier_num,
      "(minimum dissimilarity)\n"
    )
    cat(paste0(
      "Final number of outliers:\t ", backtrack_num,
      " (backtrack: max_total_rise = ", max_total_rise,
      ", max_step_rise = ", max_step_rise, ")\n"
    ))
  }

  invisible(x)
}

# ==============================================================================

#' Constructor for "outliermbc_lcwm" S3 object.
#'
#' @param x List.
#'
#' @returns "outliermbc_lcwm" S3 object.
new_outliermbc_lcwm <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class = c("outliermbc_lcwm", "list"))
}

# ------------------------------------------------------------------------------

#' Validator for `"outliermbc_lcwm"` S3 class.
#'
#' @param x List.
#'
#' @returns NULL
validate_outliermbc_lcwm <- function(x) {
  values <- unclass(x)

  obs_num <- lapply(
    values[c("labels", "outlier_bool", "outlier_rank", "gross_outs")],
    length
  )
  if (length(unique(obs_num)) != 1) {
    stop(
      paste0(
        "labels, outlier_bool, outlier_rank, and gross_outs should have the ",
        "same length."
      ),
      call. = FALSE
    )
  }

  iter_num <- c(lapply(
    values[c("loglike", "removal_dens", "distrib_diff_vec")],
    length
  ), nrow(values$distrib_diff_mat))
  if (length(unique(iter_num)) != 1) {
    stop(
      paste0(
        "loglike, removal_dens, and distrib_diff_vec should have the same ",
        "length. This should also equal the number of rows in distrib_diff_mat."
      ),
      call. = FALSE
    )
  }
}

# ------------------------------------------------------------------------------

#' plot method for `"outliermbc_lcwm"` S3 class.
#'
#' @param x List
#' @param backtrack Logical
#' @param ... Other arguments
#'
#' @returns A ggplot
#' @export
plot.outliermbc_lcwm <- function(x, backtrack = FALSE, ...) {
  plot.outliermbc_gmm(x, backtrack)
}

# ------------------------------------------------------------------------------

#' print method for `"outliermbc_lcwm"` S3 class.
#'
#' @param x List
#' @param backtrack Logical
#' @param ... Other arguments
#' @inheritParams backtrack
#'
#' @returns A ggplot
#' @export
print.outliermbc_lcwm <- function(
    x, backtrack = FALSE, max_total_rise = 0.1, max_step_rise = 0.05, ...) {
  print.outliermbc_gmm(x, backtrack, max_total_rise, max_step_rise)
}

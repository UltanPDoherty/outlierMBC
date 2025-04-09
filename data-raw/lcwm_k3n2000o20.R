lcwm_k3n2000o20 <- simulate_lcwm(
  n = c(600, 600, 800),
  mu = list(c(3), c(6), c(3)),
  sigma = list(as.matrix(1), as.matrix(0.1), as.matrix(1)),
  beta = list(c(0, 0), c(-75, 15), c(0, 5)),
  error_sd = c(1, 1, 1),
  outlier_num = c(6, 6, 8),
  outlier_type = "x_and_y",
  seed = 123,
  prob_range = c(1e-8, 1e-6),
  range_multipliers = c(1, 2)
)

usethis::use_data(lcwm_k3n2000o20, overwrite = TRUE)

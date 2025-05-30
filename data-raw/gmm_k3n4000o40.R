gmm_k3n4000o40 <- simulate_gmm(
  n = c(2000, 1000, 1000),
  mu = list(c(-1, 0), c(+1, -1), c(+1, +1)),
  sigma = list(diag(c(0.2, 4 * 0.2)), diag(c(0.2, 0.2)), diag(c(0.2, 0.2))),
  outlier_num = 40,
  seed = 123,
  crit_val = 0.9999,
  range_multiplier = 1.5
)
usethis::use_data(gmm_k3n4000o40, overwrite = TRUE)

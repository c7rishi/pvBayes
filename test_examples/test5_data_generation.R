temp <- signal_generation(
    I = 47,
    J = 7,
    sig_position = list(
      c(1,1),
      c(45,0)
    ),
    sig_strength = list(
      1,
      1.1
    ),
    sig_mean_vec = rep(3, 6),
    sig_corr_mat = matrix(0.1, nrow = 6, ncol = 6) %>% `diag<-`(1),
    var_vec = rep(.05, 6)
)

data_generation(
  lambda_mat = temp$lambda,
  signal_position = temp$signal_true,
  omega_vec = c(rep(0.1, 6),0)
)


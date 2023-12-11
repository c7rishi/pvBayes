data {

  int<lower=0> I;
  int<lower=0> J;
  array[I, J] real<lower=0> E;
  array[I, J] int<lower=0> n;

}

transformed data {

  array[I, J] real log_E = log(E);

}

parameters {

  real<lower = 0> tau;
  real<lower = 0> sigma_AE;
  vector<lower = 0>[J] sigma_Drug;

  real<lower = 0> sigma_indep2;

  cholesky_factor_corr[J] L_rho_Drug;

  vector[J] beta_Drug;

  array[I] real beta_AE;
  array[I, J] real<lower=0> theta;
  array[I, J] real log_lambda_resid;

}

transformed parameters {
/*
  cov_matrix[J] sigma_beta_Drug_joint;

  for (j1 in 1 : J) {
    for (j2 in 1 : J) {

      if (j1 == j2) {
        sigma_beta_Drug_joint[j1, j2] = sigma_Drug^2;
      } else if (j1 != J && j2 != J){
        sigma_beta_Drug_joint[j1, j2] = sigma_Drug^2 * rho_Drug[j1, j2];
      } else {
        sigma_beta_Drug_joint[j1, j2] = 0;
      }

    }
  }
*/
  vector[J] zero_mean;

  for (j in 1:J) {
    zero_mean[j] = 0;
  }

  array[I, J] real log_lambda;
  array[I, J] real log_mu;

  for (i in 1 : I){
    for (j in 1 : J ){
      log_lambda[i, j] = beta_AE[i] + beta_Drug[j] + log_lambda_resid[i, j];
      log_mu[i, j] = log_lambda[i, j] + log_E[i, j];
    }
  }




}

model {

  tau ~ cauchy(0, 1);
  //L_rho_Drug ~ uniform(-1, 1);
  sigma_AE ~ cauchy(0, 1);
  sigma_Drug ~ cauchy(0, 1);
  sigma_indep2 ~ cauchy(0, 1);
  beta_Drug ~ multi_normal_cholesky(zero_mean, diag_pre_multiply(sigma_Drug, L_rho_Drug));

  L_rho_Drug ~ lkj_corr_cholesky(1);
  //sigma ~ cauchy(0, 5);

  // target += multi_normal_lpdf(beta_Drug | zero_mean, sigma_beta_Drug_joint);

  for (i in 1 : I){
    for (j in 1 : J){
      theta[i, j] ~ cauchy (0, 1);
      log_lambda_resid[i, j] ~ normal ( 0, sqrt(tau^2 * theta[i, j]^2 +sigma_indep2^2) );
      n[i, j] ~ poisson_log ( log_mu[i, j] );
    }
    beta_AE[i] ~ normal (0, sigma_AE);
  }

}

generated quantities{

  array[I, J] int<lower=0> n_pred;
  array[I, J] real<lower=0> lambda;
  array[I, J] real<lower=0> lambda_resid;

  matrix[J,J] rho_Drug;
  rho_Drug = multiply_lower_tri_self_transpose(L_rho_Drug);

  for (i in 1 : I){
    for (j in 1 : J){
      n_pred[i, j] = poisson_log_rng ( log_mu[i, j] );
      lambda[i, j] = exp(log_lambda[i, j]);
      lambda_resid[i, j] = exp(log_lambda_resid[i, j]);
    }
  }

}

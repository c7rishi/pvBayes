// Poisson model assuming independent and specify horseshoe prior on Poisson mean
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
  real<lower = 0> sigma_Drug;
  real <lower = 0, upper = 1> rho_Drug;

  vector[J] beta_Drug;

  real<lower = 0> sigma_indep2;

  array[I] real beta_AE;
  array[I, J] real<lower=0> theta;
  array[I, J] real log_lambda_resid;

}

transformed parameters {

  cov_matrix[J] sigma_beta_Drug_joint;

  for (j1 in 1 : J) {
    for (j2 in 1 : J) {

      if (j1 == j2) {
        sigma_beta_Drug_joint[j1, j2] = sigma_Drug^2;
      } else if (j1 != J && j2 != J){
        sigma_beta_Drug_joint[j1, j2] = sigma_Drug^2 * rho_Drug;
      } else {
        sigma_beta_Drug_joint[j1, j2] = 0;
      }

    }
  }

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
  rho_Drug ~ uniform(-1, 1);
  sigma_AE ~ cauchy(0, 1);
  sigma_Drug ~ cauchy(0, 1);
  sigma_indep2 ~ cauchy(0,1);
  target += multi_normal_lpdf(beta_Drug | zero_mean, sigma_beta_Drug_joint);

  for (i in 1 : I){
    for (j in 1 : J){
      theta[i, j] ~ cauchy (0, 1);
      log_lambda_resid[i, j] ~ normal ( 0, sqrt(sigma_indep2^2 + tau^2 * theta[i, j]^2) );
      n[i, j] ~ poisson_log ( log_mu[i, j] );
    }
    beta_AE[i] ~ normal (0, sigma_AE);
  }

}

generated quantities{

  array[I, J] int<lower=0> n_pred;
  array[I, J] real<lower=0> lambda;
  array[I, J] real<lower=0> lambda_resid;
  for (i in 1 : I){
    for (j in 1 : J){
      n_pred[i, j] = poisson_log_rng ( log_mu[i, j] );
      lambda[i, j] = exp(log_lambda[i, j]);
      lambda_resid[i, j] = exp(log_lambda_resid[i, j]);
    }
  }

}

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
  real <lower = 0, upper = 1> rho_AE;

  vector[I] beta_AE;

  array[J] real beta_Drug;
  array[I, J] real<lower=0> theta;
  array[I, J] real log_lambda_resid;

}

transformed parameters {

  cov_matrix[I] sigma_beta_AE_joint;

  for (i1 in 1 : I) {
    for (i2 in 1 : I) {

      if (i1 == i2) {
        sigma_beta_AE_joint[i1, i2] = sigma_AE^2;
      } else if (i1 != I && i2 != I){
        sigma_beta_AE_joint[i1, i2] = sigma_AE^2 * rho_AE;
      } else {
        sigma_beta_AE_joint[i1, i2] = 0;
      }

    }
  }

  vector[I] zero_mean;

  for (i in 1:I) {
    zero_mean[i] = 0;
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
  rho_AE ~ uniform(-1, 1);
  sigma_AE ~ cauchy(0, 1);
  sigma_Drug ~ cauchy(0, 1);

  target += multi_normal_lpdf(beta_AE | zero_mean, sigma_beta_AE_joint);


  for (j in 1 : J){
    for (i in 1 : I){
      theta[i, j] ~ cauchy (0, 1);
      log_lambda_resid[i, j] ~ normal ( 0, tau * theta[i, j] );
      n[i, j] ~ poisson_log ( log_mu[i, j] );
    }
    beta_Drug[j] ~ normal (0, sigma_Drug);
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

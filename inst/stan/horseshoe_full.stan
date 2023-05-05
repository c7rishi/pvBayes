data {

  int<lower=0> I;
  int<lower=0> J;
  array[I, J] real<lower=0> E;
  array[I, J] int<lower=0> n;

}

transformed data {

  array[I, J] real log_E;
  log_E = log(E);
}

parameters {

  real<lower=0> tau;
  real log_beta_0;
  array[I] real log_beta_AE;
  array[J] real log_beta_Drug;
  real<lower=0> sigma_beta_0;
  real<lower=0> sigma_beta_AE;
  real<lower=0> sigma_beta_Drug;
  array[I, J] real<lower=0> log_lambda;
  array[I, J] real<lower=0> theta;

}

transformed parameters{

  array[I, J] real<lower=0> lambda;

  array[I, J] real log_mu;

  for (i in 1 : I){
    for (j in 1 :J ){
      log_mu[i, j] = log_lambda[i, j] + log_beta_0 + log_beta_AE[i] + log_beta_Drug[j] + log_E[i, j];
      lambda[i, j] = exp(log_lambda[i, j]);
    }
  }


}

model {

  tau ~ cauchy(0, 1);

  sigma_beta_0 ~ cauchy(0, 1);
  sigma_beta_AE ~ cauchy(0, 1);
  sigma_beta_Drug ~ cauchy(0, 1);

  log_beta_0 ~ normal(0, sigma_beta_0);
  log_beta_AE ~ normal(0, sigma_beta_AE);
  log_beta_Drug ~ normal(0, sigma_beta_Drug);

  for (i in 1 : I){
    for (j in 1 : J){
      theta[i,j] ~ cauchy (0, 1);
      log_lambda[i, j] ~ normal(0, tau * theta[i, j]);

      n[i, j] ~ poisson_log( log_mu[i, j] );
    }
  }

}

generated quantities{

  array[I, J] int<lower=0> n_pred;

  for (i in 1 : I){
    for (j in 1 : J){
      n_pred[i, j] = poisson_log_rng ( log_mu[i, j] );
    }
  }

}




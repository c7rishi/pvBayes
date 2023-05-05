data {

  int<lower=0> I;
  int<lower=0> J;
  array[I, J] real<lower=0> E;
  array[I, J] int<lower=0> n;

}


parameters {

  real<lower=0> tau;
  real<lower = 0> sigma_AE;
  real<lower = 0> sigma_Drug;
  real beta_0;
  array[I] real beta_AE_raw_unscaled;
  array[J] real beta_Drug_raw_unscaled;
  array[I, J] real<lower=0> log_lambda;
  array[I, J] real<lower=0> theta;

}

transformed parameters{

  real mean_beta_AE = mean(beta_AE_raw_unscaled);
  real mean_beta_Drug = mean(beta_Drug_raw_unscaled);

  array[I] real beta_AE;
  array[J] real beta_Drug;
  array[I, J] real log_mu;

  for (i in 1: I){
    beta_AE[i] = sigma_AE * (beta_AE_raw_unscaled[i] - mean_beta_AE);
  }

  for (j in 1 : J) {
    beta_Drug[j] = sigma_Drug * (beta_Drug_raw_unscaled[j] - mean_beta_Drug);
  }

  for (i in 1 : I){
    for (j in 1 : J){
      log_mu[i, j] = beta_0 + beta_AE[i] + beta_Drug[j] + log_lambda[i, j];
    }
  }


}

model {

  tau ~ cauchy(0, 1);
  beta_0 ~ normal(0, 10);
  beta_AE ~ normal(0, 10);
  beta_Drug ~ normal(0, 10);
  sigma_AE ~ uniform(1, 2);
  sigma_Drug ~ uniform(1, 2);

  for (i in 1 : I){
    for (j in 1 : J){
      theta[i, j] ~ cauchy (0, 1);
      log_lambda[i, j] ~ normal ( 0, tau * theta[i, j] );
      n[i, j] ~ poisson_log( log_mu[i, j] );
    }
  }

}

generated quantities{

  array[I, J] int<lower=0> n_pred;
  array[I, J] real<lower=0> lambda;

  for (i in 1 : I){
    for (j in 1 : J){
      n_pred[i, j] = poisson_log_rng ( log_mu[i, j] );
      lambda[i, j] = exp(log_lambda[i, j]);
    }
  }

}




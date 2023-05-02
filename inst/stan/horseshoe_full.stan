data {

  int<lower=0> I;
  int<lower=0> J;
  array[I, J] real<lower=0> E;
  array[I, J] int<lower=0> n;

}


parameters {

  array[I, J] real<lower=0> lambda;
  real<lower=0> tau;
  array[I, J] real<lower=0> theta;

  array[I] real beta_AE_raw_unscaled;
  array[J] real beta_Drug_raw_unscaled;
  real beta_0;

  real<lower = 0> sigma_AE;
  real<lower = 0> sigma_Drug;


}

transformed parameters{

  array[I] real beta_AE;
  array[J] real beta_Drug;
  array[I, J] real log_mu;
  array[I, J] real log_lambda;

  real mean_beta_AE = sum(beta_AE_raw_unscaled)/I;
  real mean_beta_Drug = sum(beta_Drug_raw_unscaled)/J;

  log_lambda = log(lambda);
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



  beta_0 ~ normal(0, 10);
  beta_AE ~ normal(0, 10);
  sigma_AE ~ uniform(1, 2);
  beta_Drug ~ normal(0, 10);
  sigma_Drug ~ uniform(1, 2);

  for (i in 1 : I){
    for (j in 1 : J){
      n[i, j] ~ poisson_log(log_mu[i, j]);
    }
  }

  tau ~ cauchy(0, 1);

  for (i in 1 : I){
    for(j in 1 : J){
      theta[i, j] ~ cauchy (0, 1);
      log_lambda[i, j] ~ normal ( 0, tau * theta[i, j] );
    }
  }


}

generated quantities{

  array[I, J] int<lower=0> n_pred;
  for (i in 1 : I){
    for (j in 1 : J){
      n_pred[i, j] = poisson_log_rng(log_mu[i, j]);
    }
  }

}



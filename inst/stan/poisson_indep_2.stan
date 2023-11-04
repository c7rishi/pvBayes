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
  real<lower = 0> sigma_indep;

  array[I, J] real<lower=0> theta;
  array[I, J] real log_lambda_indep;
  array[I, J] real log_lambda_resid;

}

transformed parameters {

  array[I, J] real log_lambda;
  array[I, J] real log_mu;

  for (i in 1 : I){
    for (j in 1 :J ){
      log_lambda[i, j] = log_lambda_indep[i, j] + log_lambda_resid[i, j];
      log_mu[i, j] = log_lambda[i, j] + log_E[i, j];
      //log_lambda[i,j] = log_mu[i,j] - log_E[i,j];
    }
  }

}

model {

  tau ~ cauchy(0, 1);
  sigma_indep ~ cauchy(0, 1);

  for (i in 1 : I){
    for (j in 1 : J){
      theta[i, j] ~ cauchy (0, 1);
      log_lambda_indep[i, j] ~ normal ( 0, sqrt(sigma_indep^2) );
      log_lambda_resid[i, j] ~ normal ( 0, sqrt(tau^2 * theta[i, j]^2));
      n[i, j] ~ poisson_log ( log_mu[i, j] );
    }
  }

}

generated quantities{

  array[I, J] int<lower=0> n_pred;
  array[I, J] real<lower=0> lambda;
  array[I, J] real<lower=0> lambda_indep;
  array[I, J] real<lower=0> lambda_resid;
  for (i in 1 : I){
    for (j in 1 : J){
      n_pred[i, j] = poisson_log_rng ( log_mu[i, j] );
      lambda[i, j] = exp(log_lambda[i, j]);
      lambda_indep[i, j] = exp(log_lambda_indep[i, j]);
      lambda_resid[i, j] = exp(log_lambda_resid[i, j]);
    }
  }

}

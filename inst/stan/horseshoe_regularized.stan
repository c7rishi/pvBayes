// Poisson model with regularized horseshoe prior ( Piironen 2017)
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

  real<lower=0> tau;
  real<lower=0> c2;

  array[I, J] real<lower=0> log_lambda;
  array[I, J] real<lower=0> theta;

}

transformed parameters {



}


model {

  array[I, J] real theta_tilde;

  c2 ~ inv_gamma(2, 8);

  tau ~ cauchy(0, 1);

  for (i in 1 : I){
    for (j in 1 : J){
      theta[i, j] ~ cauchy (0, 1);
      theta_tilde[i, j] = c2 * theta[i, j]^2 / (c2 + tau^2 * theta[i, j]^2);
      log_lambda[i, j] ~ normal ( 0, tau * theta_tilde[i, j] );
      n[i, j] ~ poisson_log ( log_lambda[i, j] + log_E[i, j] );
    }
  }

}

generated quantities{

  array[I, J] int<lower=0> n_pred;
  array[I, J] real<lower=0> lambda;

  for (i in 1 : I){
    for (j in 1 : J){
      n_pred[i, j] = poisson_log_rng ( log_lambda[i, j] + log_E[i, j] );
      lambda[i, j] = exp(log_lambda[i, j]);
    }
  }

}



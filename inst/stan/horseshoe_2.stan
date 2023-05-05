// Horseshoe prior with restriction on local means
data {

  int<lower=0> I;
  int<lower=0> J;
  array[I, J] real<lower=0> E;
  array[I, J] int<lower=0> n;

  int<lower=0> gamma;

}

transformed data {

  array[I, J] real log_E_ratio;
  array[I, J] real log_E = log(E);

  for(i in 1 : I){
    for(j in 1 : J){
      log_E_ratio[i, j] = log(E[i, j])-log(1+E[i, j]);
    }

  }

}

parameters {

  real<lower=0> tau;
  array[I, J] real log_lambda;
  array[I, J] real<lower=0> theta;

}

model {

  tau ~ cauchy(0, 1);

  for (i in 1 : I){
    for (j in 1 : J){
      theta[i, j] ~ cauchy (0, exp(gamma * log_E_ratio[i, j]));
      log_lambda[i, j] ~ normal ( 0, tau * theta[i, j] );
      n[i, j] ~ poisson_log( log_lambda[i, j] + log_E[i, j]);
    }
  }

}

generated quantities {

  array[I, J] int<lower=0> n_pred;
  array[I, J] real<lower=0> lambda;

  for (i in 1 : I){
    for (j in 1 : J){
      n_pred[i, j] = poisson_log_rng ( log_lambda[i, j] + log_E[i, j] );
      lambda[i, j] = exp(log_lambda[i, j]);
    }
  }

}

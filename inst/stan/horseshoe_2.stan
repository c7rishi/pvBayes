data {

  int<lower=0> I;
  int<lower=0> J;
  array[I, J] real<lower=0> E;
  array[I, J] int<lower=0> n;
  real<lower=0> b;

}

transformed data {

  array[I, J] real log_E_ratio;
  array[I, J] real log_E;

  for(i in 1 : I){
    for(j in 1 : J){
      log_E_ratio[i, j] = log(E[i, j])-log(1+E[i, j]);
    }

  }

  log_E = log(E);

}

parameters {

  array[I, J] real<lower=0> lambda;
  real<lower=0> tau;
  array[I, J] real<lower=0> theta;
  real<lower=0> gamma;
}

model {

  array[I, J] real log_lambda;
  log_lambda = log(lambda);

  gamma ~ uniform(0, b);

  for (i in 1 : I){
    for (j in 1 : J){
      n[i, j] ~ poisson_log( log_lambda[i, j] + log_E[i, j]);
    }
  }

  tau ~ cauchy(0, 1);

  for (i in 1 : I){
    for(j in 1 : J){
      theta[i, j] ~ cauchy (0, exp(gamma * log_E_ratio[i, j]));
      log_lambda[i, j] ~ normal ( 0, tau * theta[i, j] );
    }
  }

}


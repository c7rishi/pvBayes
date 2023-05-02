data {

  int<lower=0> I;
  int<lower=0> J;
  array[I, J] real<lower=0> E;
  array[I, J] int<lower=0> n;
  //real<lower=0> gamma;
}

/*transformed data {

  array[I, J] real E_ratio;

  for(i in 1 : I){
    for(j in 1 : J){
      E_ratio[i, j] = E[i, j]/(1+E[i, j]);
    }
  }

}*/

parameters {

  array[I, J] real<lower=0> lambda;
  real<lower=0> tau;
  array[I, J] real<lower=0> theta;

}

model {

  array[I, J] real log_lambda;
  log_lambda = log(lambda);

  for (i in 1 : I){
    for (j in 1 : J){
      n[i, j] ~ poisson ( lambda[i, j] * E[i, j] );
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
      n_pred[i, j] = poisson_rng ( lambda[i, j] * E[i, j] );
    }
  }
}

data {

  int<lower=0> I;
  int<lower=0> J;
  array[I, J] real<lower=0> E;
  array[I, J] int<lower=0> n;
  //real<lower=0> gamma;
}

transformed data {

  array[I, J] real log_E;

  for(i in 1 : I){
    for(j in 1 : J){
      log_E[i, j] = log(E[i, j]);
    }
  }

}

parameters {

  real<lower=0> tau;
  array[I, J] real<lower=0> theta;
  real<lower = 0> sigma_indep;
  array[I, J] real log_mu;

}

transformed parameters {

  array[I, J] real log_lambda;
  array[I, J] real<lower=0> lambda;

  for (i in 1 : I){
    for (j in 1 :J ){
      log_lambda[i,j] = log_mu[i,j] - log_E[i,j];
      lambda[i,j] = exp(log_lambda[i,j]);
    }
  }



}

model {


  for (i in 1 : I){
    for (j in 1 : J){
      n[i, j] ~ poisson_log ( log_mu[i, j] );
    }
  }

  tau ~ cauchy(0, 1);
  sigma_indep ~ cauchy(0, 1);
  for (i in 1 : I){
    for(j in 1 : J){
      theta[i, j] ~ cauchy (0, 1);
      log_mu[i, j] ~ normal ( log_E[i, j], sqrt(sigma_indep^2+tau^2 * theta[i, j]^2) );
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

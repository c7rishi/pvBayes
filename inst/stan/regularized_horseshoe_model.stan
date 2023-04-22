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

  array[I, J] real<lower=0> lambda;
  real<lower=0> tau;
  array[I, J] real<lower=0> theta;
  real<lower=0> c;
}

transformed parameters {

  real<lower=0> c2;
  c2 = c^2;


}


model {

  array[I, J] real log_lambda;
  log_lambda = log(lambda);



  array[I, J] real theta_tilde;


  for (i in 1 : I){
    for (j in 1 : J){
      n[i, j] ~ poisson ( lambda[i, j] * E[i, j] );
    }
  }

  c2 ~ inv_gamma(2, 8);

  tau ~ cauchy(0, 1);

  for (i in 1 : I){
    for(j in 1 : J){
      theta_tilde[i, j] = c2 * theta[i, j]^2 / (c2 + tau^2 * theta[i, j]^2);
      theta[i, j] ~ cauchy (0, 1);
      log_lambda[i, j] ~ normal ( 0, tau * theta_tilde[i, j] );
    }
  }

}


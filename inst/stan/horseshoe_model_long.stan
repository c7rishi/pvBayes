data {
  
  int<lower=0> N;
  array[N] real<lower=0> E;
  array[N] int<lower=0> n;

}

transformed data {
  
  array[N] real log_E;
  log_E = log(E);

}

parameters {
  
  array[N] real<lower=0> lambda;
  real<lower=0> tau;
  array[N] real<lower=0> theta;

}

transformed parameters {
  
  array[N] real log_lambda;
  log_lambda = log(lambda);

}

model {
  
  for (i in 1 : N){
    
    n[i] ~ poisson ( lambda[i] * E[i] );
  
  }
  
  tau ~ cauchy(0, 1);
  
  for (i in 1 : N){
    
    theta[i] ~ cauchy (0, 1);
    log_lambda[i] ~ normal ( 0, tau * theta[i] );
  
  }
  
}


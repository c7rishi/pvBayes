functions {
  
  real lhgf_2F1(real a, real b, real c, real z){
    
    real tau;
    real y_hat;
    real r_2_1;
    
    tau = z * (a - b) - c;
    
    y_hat = 2 * b / (sqrt( tau^2 - 4 * b * z * (c - a) ) - tau);
    
    r_2_1 = y_hat^2 / b + (1 - y_hat)^2 / (c - b) 
    ./ - a * z^2 / (1 - z*y_hat)^2 * y_hat^2 / b * (1 - y_hat)^2 / (c - b);
    
    return (c - 0.5)*log(c) - 0.5*log(r_2_1) + b*log(y_hat) - b*log(b) 
    ./ + (c - b)*log(1 - y_hat) - (c - b)*log(c - b) - a*log(1 - z * y_hat);
    
  } 
  
  real GH_lpdf(real kappa, real a, real b, real phi, real gamma){
    
    real out;
    
    out = - lbeta(a, b) - lhgf_2F1(gamma, a, a+b, 1-phi)
    ./ + (a-1) * log(kappa) + (b-1) * log(1-kappa)
    ./ - gamma * log(1-(1-phi)*kappa);
    
    return out;
    
  }
  
}

data {
  
  int<lower=0> N;
  array[N] real<lower=0> E;
  array[N] int<lower=0> n;
  
}

parameters {
  
  array[N] real<lower=0, upper=1> kappa;
  real<lower=0> alpha;
  real<lower=0> tau;

}

transformed parameters{
  
  array[N] real<lower=0, upper=1> kappa_tilde;
  array[N] real<lower=0> beta;
  
  for (i in 1:N){
    
    kappa_tilde[i] = 1 / (1 + E[i] * (1 / kappa[i]-1));
    beta[i] =  exp(log( kappa_tilde[i] ) - log( 1 - kappa_tilde[i] ));
    
  }

}

model {
  
  for (i in 1 : N){
    
    n[i] ~ neg_binomial ( alpha, beta[i] );
    
  }
  
  alpha ~ cauchy(0,10);
  tau ~ cauchy(0,10);
  
  for (i in 1 : N){
    
    target += GH_lpdf(kappa[i] | 0.5, 0.5, tau*tau, 0.5);
    
  }
}

generated quantities {
  
  array[N] real<lower=0> lambda;
  
  for (i in 1 : N){
    lambda[i] = gamma_rng(n[i] + alpha, E[i] / (1-kappa_tilde[i]));
  }
  
}


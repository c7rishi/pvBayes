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
  array[N] real<lower=0, upper=1> omega;
}

transformed parameters {

  array[N] real log_lambda;
  log_lambda = log(lambda);

}

model {

  for (i in 1 : N) {
    if (n[i] == 0) {
      target += log_sum_exp(bernoulli_lpmf(1 | omega[i]),
                            bernoulli_lpmf(0 | omega[i])
                              + poisson_lpmf(0 | lambda[i]*E[i]));
    } else {
      target += bernoulli_lpmf(0 | omega[i])
                  + poisson_lpmf(n[i] | lambda[i]*E[i]);
    }
  }

  tau ~ cauchy(0, 1);

  for (i in 1 : N){

    omega[i] ~ beta(0.5, 0.5);
    theta[i] ~ cauchy (0, 1);
    log_lambda[i] ~ normal ( 0, tau * theta[i] );

  }

}


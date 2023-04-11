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
  array[J] real<lower=0, upper=1> omega;

}

transformed parameters {

  array[I, J] real log_lambda;
  log_lambda = log(lambda);

}

model {

  for (i in 1 : I) {
    for (j in 1 : J) {
      if (n[i,j] == 0) {

        target += log_sum_exp(bernoulli_lpmf(1 | omega[j]),
                               bernoulli_lpmf(0 | omega[j])
                                + poisson_lpmf(0 | lambda[i,j]*E[i,j]));

      } else {

        target += bernoulli_lpmf(0 | omega[j])
                   + poisson_lpmf(n[i,j] | lambda[i,j]*E[i,j]);
      }

    }
  }

  tau ~ cauchy(0, 1);

  for (j in 1 : J){

    omega[j] ~ beta(0.5, 0.5); #Jeffrey's prior

    for (i in 1 : I){

      theta[i,j] ~ cauchy (0, 1);
      log_lambda[i,j] ~ normal ( 0, tau * theta[i,j] );

    }
  }
}


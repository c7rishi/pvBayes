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

  real<lower=0> tau;
  array[J] real<lower=0, upper=1> omega;
  array[I, J] real log_lambda_raw;
  array[I, J] real<lower=0> theta;
  real<lower = 0> sigma_indep2;

}

transformed parameters {

  array[I, J] real<lower=0> lambda_raw;

  for (j in 1 : J){
    for (i in 1 : I){
      lambda_raw[i, j] = exp(log_lambda_raw[i, j]);
    }
  }


}

model {

  tau ~ cauchy(0, 1);
  omega ~ beta(0.5, 0.5); //Jeffrey's prior

  for (i in 1 : I) {
    for (j in 1 : J) {

      theta[i,j] ~ cauchy (0, 1);
      log_lambda_raw[i,j] ~ normal ( 0, sqrt(tau^2 * theta[i, j]^2 +sigma_indep2^2) );

      if (n[i, j] == 0) {

        target += log_sum_exp(bernoulli_lpmf(1 | omega[j]),
        bernoulli_lpmf(0 | omega[j])
        + poisson_lpmf(0 | lambda_raw[i, j] * E[i, j] ) );

      } else {

        target += bernoulli_lpmf(0 | omega[j])
        + poisson_lpmf(n[i, j] | lambda_raw[i, j] * E[i, j] );

      }

    }
  }

}

generated quantities {

  array[I, J] real<lower=0, upper=1> zi;
  array[I, J] int<lower=0> n_pred;
  array[I, J] int<lower=0, upper=1> zi_pred;
  array[I, J] real<lower=0> lambda;

  for (j in 1 : J){
    for (i in 1 : I){

      n_pred[i, j] = poisson_log_rng ( log_lambda_raw[i, j] + log_E[i, j] );

      if (n[i, j] == 0) {

        zi[i, j] = omega[j] / ( omega[j] + (1-omega[j]) * exp(-lambda_raw[i, j] * E[i, j]) );

      } else{

        zi[i, j] = 0;

      }
      zi_pred[i, j] = bernoulli_rng( zi[i, j] );
      lambda[i, j] = (1 - zi_pred[i, j]) * lambda_raw[i, j];
    }
  }

}

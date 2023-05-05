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
  real<lower = 0> sigma_indep;
  array[J] real<lower=0, upper=1> omega;
  array[I, J] real log_mu;
  array[I, J] real<lower=0> theta;

}

transformed parameters {

  array[I, J] real<lower=0> lambda;
  array[I, J] real log_lambda;

  for (i in 1 : I){
    for (j in 1 :J ){
      log_lambda[i,j] = log_mu[i,j] - log_E[i,j];
      lambda[i, j] = exp(log_lambda[i, j]);
    }
  }

}

model {

  tau ~ cauchy(0, 1);
  omega ~ beta(0.5, 0.5); //Jeffrey's prior
  sigma_indep ~ cauchy(0, 1);

  for (i in 1 : I) {
    for (j in 1 : J) {

      theta[i,j] ~ cauchy (0, 1);
      log_mu[i, j] ~ normal ( log_E[i, j], sqrt(sigma_indep^2 + tau^2 * theta[i, j]^2) );

      if (n[i, j] == 0) {

        target += log_sum_exp(bernoulli_lpmf(1 | omega[j]),
        bernoulli_lpmf(0 | omega[j])
        + poisson_log_lpmf(0 | log_mu[i, j] ) );

      } else {

        target += bernoulli_lpmf(0 | omega[j])
        + poisson_log_lpmf(n[i, j] | log_mu[i, j] );

      }

    }
  }

}

generated quantities {

  array[I, J] real<lower=0, upper=1> zi;
  array[I, J] int<lower=0> n_pred;

  for (j in 1 : J){
    for (i in 1 : I){

      n_pred[i, j] = poisson_log_rng ( log_mu[i, j] );

      if (n[i, j] == 0) {

        zi[i, j] = omega[j] / ( omega[j] + (1-omega[j]) * exp(-lambda[i, j] * E[i, j]) );

      } else{

        zi[i, j] = 0;

      }
    }
  }

}

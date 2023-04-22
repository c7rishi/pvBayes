data {

  int<lower=0> I;
  int<lower=0> J;
  array[I, J] real<lower=0> E;
  array[I, J] int<lower=0> n;
  int<lower=1> t_nu;
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
  real<lower=0> c;

}

transformed parameters {

  array[I, J] real log_lambda;
  log_lambda = log(lambda);

  real<lower=0> c2;
  c2 = c^2;

}

model {

  array[I, J] real theta_tilde;

  c2 ~ inv_gamma(2, 8);

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
  omega ~ beta(0.5, 0.5); //Jeffrey's prior

  for (j in 1 : J){
    for (i in 1 : I){
      theta[i, j] ~ student_t (t_nu, 0, 1);
      //theta_tilde[i, j] = c2 * theta[i, j]^2 / (c2 + tau^2 * theta[i, j]^2);
      //theta[i,j] ~ cauchy (0, 1);
      //log_lambda[i, j] ~ normal ( 0, tau * theta_tilde[i, j] );
      log_lambda[i, j] ~ normal ( 0, tau* theta[i, j] );

      //
    }
  }
}

generated quantities {

  array[I, J] real<lower=0, upper=1> zi;

  for (j in 1 : J){
    for (i in 1 : I){
      if (n[i, j] == 0) {

        zi[i, j] = omega[j] / ( omega[j] + (1-omega[j])*exp(-lambda[i, j]*E[i, j]) );

      } else{

        zi[i, j] = 0;

      }
    }
  }

}

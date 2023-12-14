data {

  int<lower=0> I;
  int<lower=0> J;
  array[I, J] real<lower=0> E;
  array[I, J] int<lower=0> n;

}

transformed data {

  array[I, J] real log_E = log(E);

  vector[J-1] zero_mean;

  for (j in 1:(J-1)) {
    zero_mean[j] = 0;
  }

}

parameters {


  //array[I, J] real beta_indep;

  //real<lower = 0> tau_indep;
  //array[I, J] real<lower=0> theta_indep;
  //real<lower = 0> sigma_indep2_indep;

  real<lower = 0> tau;
  matrix<lower=0>[I, J-1] theta;
  real<lower = 0> sigma_indep2;

  cholesky_factor_corr[J-1] L_rho_Drug;


  matrix<lower=0>[I, J-1] log_lambda_relevant;
  array[I] real<lower=0> log_lambda_other;

  array[J] real<lower=0, upper=1> omega;


}

transformed parameters {


  array[I, J] real log_lambda;
  array[I, J] real log_mu;

  matrix<lower=0>[I, J-1] sigma_lambda_relevant;
  //array[I] real<lower=0> sigma_lambda_other;

  for (i in 1 : I) {

    for (j in 1 : (J-1) ) {

      log_lambda[i, j] = log_lambda_relevant[i, j];
      sigma_lambda_relevant[i, j] = sqrt(tau^2 * theta[i, j]^2 + sigma_indep2^2);
      //sigma_lambda[i, j] = sigma_lambda_relevant[i, j];

    }

    log_lambda[i, J] = log_lambda_other[i];
    //sigma_lambda_other[i] = sqrt(tau^2 * theta[i, J]^2 + sigma_indep2^2);
    //sigma_lambda[i, J] = sigma_lambda_other[i];

  }


  for (i in 1 : I){

    for (j in 1 : J ){

      //log_mu[i, j] = beta_indep[i, j] + log_lambda[i, j] + log_E[i,j];
      log_mu[i, j] = log_lambda[i, j] + log_E[i, j];

    }
  }


}

model {

  //tau_indep ~ cauchy(0, 1);
  //sigma_indep2_indep ~ cauchy(0, 1);

  tau ~ cauchy(0, 1);
  sigma_indep2 ~ cauchy(0, 1);

  L_rho_Drug ~ lkj_corr_cholesky(1);

  for (i in 1 : I){

    log_lambda_relevant[i] ~ multi_normal_cholesky(zero_mean, diag_pre_multiply(sigma_lambda_relevant[i], L_rho_Drug));
    log_lambda_other[i] ~ normal (0, sigma_indep2);

    theta[i] ~ cauchy (0, 1);

    for (j in 1 : J){

      //theta_indep[i, j] ~ cauchy (0, 1);
      //beta_indep[i, j] ~ normal (0, sqrt(tau_indep^2 * theta_indep[i, j]^2 + sigma_indep2_indep^2) );



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
  array[I, J] int<lower=0, upper=1> zi_pred;
  array[I, J] real<lower=0> lambda;

  matrix[J-1,J-1] rho_Drug;
  rho_Drug = multiply_lower_tri_self_transpose(L_rho_Drug);


  for (j in 1 : J){
    for (i in 1 : I){

      n_pred[i, j] = poisson_log_rng ( log_mu[i, j] );

      if (n[i, j] == 0) {

        zi[i, j] = omega[j] / ( omega[j] + (1-omega[j]) * exp(-exp(log_mu[i, j])) );

      } else{

        zi[i, j] = 0;

      }
      zi_pred[i, j] = bernoulli_rng( zi[i, j] );
      lambda[i, j] = (1 - zi_pred[i, j]) * exp(log_lambda[i, j]);
    }
  }

}

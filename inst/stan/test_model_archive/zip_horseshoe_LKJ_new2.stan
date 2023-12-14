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

  real beta_global;
  array[I] real beta_AE;
  array[J] real beta_Drug;

  //real<lower = 0> sigma_AE;
  //real<lower = 0> sigma_Drug;

  real<lower = 0> tau_AE;
  real<lower = 0> tau_Drug;

  real<lower = 0> sigma_indep2_AE;
  real<lower = 0> sigma_indep2_Drug;

  array[I] real<lower=0> theta_AE;
  array[J] real<lower=0> theta_Drug;

  real<lower = 0> tau;
  //vector<lower = 0>[J-1] sigma_Drug;
  array[I, J] real<lower=0> theta;
  real<lower = 0> sigma_indep2;

  cholesky_factor_corr[J-1] L_rho_Drug;


  //vector[J-1] beta_Drug_relevant;
  //real beta_Drug_other;

  matrix<lower=0>[I, J-1] log_lambda_relevant;
  //array[I, J-1] real<lower=0> log_lambda_relevant;
  array[I] real<lower=0> log_lambda_other;

  array[J] real<lower=0, upper=1> omega;

  //array[I, J] real log_lambda;

}

transformed parameters {


  array[I, J] real log_lambda;
  array[I, J] real log_mu;
  array[I, J] real log_lambda_final;
  //array[I, J] real sigma_lambda;

  matrix<lower=0>[I, J-1] sigma_lambda_relevant;
  //array[I, J-1] real<lower=0> sigma_lambda_relevant;
  array[I] real<lower=0> sigma_lambda_other;

  for (i in 1 : I) {

    for (j in 1 : (J-1) ) {

      log_lambda[i, j] = log_lambda_relevant[i, j];
      sigma_lambda_relevant[i, j] = sqrt(tau^2 * theta[i, j]^2 + sigma_indep2^2);
      //sigma_lambda[i, j] = sigma_lambda_relevant[i, j];

    }

    log_lambda[i, J] = log_lambda_other[i];
    sigma_lambda_other[i] = sqrt(tau^2 * theta[i, J]^2 + sigma_indep2^2);
    //sigma_lambda[i, J] = sigma_lambda_other[i];

  }


  for (i in 1 : I){

    for (j in 1 : J ){

      log_mu[i, j] = beta_global + beta_AE[i] + beta_Drug[j] + log_lambda[i, j] + log_E[i,j];

      // sigma_lambda[i, j] = sqrt(tau^2 * theta[i, j]^2 +sigma_indep2^2);

      log_lambda_final[i, j] = log_mu[i, j] - log_lambda[i, j];

    }
  }


}

model {

  //sigma_AE ~ cauchy(0, 1);
  //sigma_Drug ~ cauchy(0, 1);

  tau_AE ~ cauchy(0, 1);
  tau_Drug ~ cauchy(0, 1);

  sigma_indep2_AE ~ cauchy(0, 1);
  sigma_indep2_Drug ~ cauchy(0, 1);

  beta_global ~ normal(0,10);
  for ( i in 1 : I ) {
    theta_AE[i] ~ cauchy (0, 1);
    beta_AE[i] ~ normal (0, sqrt(tau_AE^2 * theta_AE[i]^2 + sigma_indep2_AE^2) );
  }
  for ( j in 1 : J ) {
    theta_Drug[j] ~ cauchy (0, 1);
    beta_Drug[j] ~ normal(0, sqrt(tau_Drug^2 * theta_Drug[j]^2 + sigma_indep2_Drug^2) );
  }


  tau ~ cauchy(0, 1);
  sigma_indep2 ~ cauchy(0, 1);

  //beta_Drug_relevant ~ multi_normal_cholesky(zero_mean, diag_pre_multiply(sigma_Drug, L_rho_Drug));

  L_rho_Drug ~ lkj_corr_cholesky(1);
  //sigma_lambda ~ multi_normal_cholesky(zero_mean, diag_pre_multiply(sigma_lambda, L_rho_Drug));


  for (i in 1 : I){

    log_lambda_relevant[i] ~ multi_normal_cholesky(zero_mean, diag_pre_multiply(sigma_lambda_relevant[i], L_rho_Drug));

    log_lambda_other[i] ~ normal ( 0, sigma_lambda_other[i] );

    for (j in 1 : J){

      theta[i, j] ~ cauchy (0, 1);
      //log_lambda[i,j] ~ normal ( 0, sqrt(tau^2 * theta[i, j]^2 +sigma_indep2^2) );

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
      lambda[i, j] = (1 - zi_pred[i, j]) * exp(log_lambda_final[i, j]);
    }
  }

}

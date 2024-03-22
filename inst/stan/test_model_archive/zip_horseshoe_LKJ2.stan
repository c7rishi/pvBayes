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

  real<lower = 0> tau;
  real<lower = 0> sigma_AE;
  vector<lower = 0>[J-1] sigma_Drug;


  //real<lower = 0> sigma_ridge;

  cholesky_factor_corr[J-1] L_rho_Drug;


  vector[J-1] beta_Drug_relevant;
  real beta_Drug_other;



  array[I] real beta_AE;
  array[I, J] real<lower=0> theta;

  array[J] real<lower=0, upper=1> omega;
  array[I, J] real log_lambda_tilde;
}

transformed parameters {


  array[I, J] real log_lambda_tilde_total;
  array[I, J] real log_mu;

  vector[J] beta_Drug;

  for (j in 1 : (J-1)) {
    beta_Drug[j] = beta_Drug_relevant[j];
  }

  beta_Drug[J] = beta_Drug_other;

  for (i in 1 : I){

    for (j in 1 : J ){
      log_lambda_tilde_total[i, j] = beta_AE[i] + beta_Drug[j] + log_lambda_tilde[i, j];
      log_mu[i, j] = log_lambda_tilde_total[i, j] + log_E[i, j];
    }
  }




}

model {

  tau ~ cauchy(0, 1);
  sigma_AE ~ cauchy(0, 1);
  sigma_Drug ~ cauchy(0, 1);
  //sigma_ridge ~ cauchy(0, 1);
  beta_Drug_other ~ normal(0, 10);
  beta_Drug_relevant ~ multi_normal_cholesky(zero_mean, diag_pre_multiply(sigma_Drug, L_rho_Drug));

  L_rho_Drug ~ lkj_corr_cholesky(1);


  for (i in 1 : I){

    for (j in 1 : J){

      theta[i, j] ~ cauchy (0, 1);
      //log_lambda_tilde[i,j] ~ normal ( 0, sqrt(tau^2 * theta[i, j]^2 +sigma_ridge^2) );
      log_lambda_tilde[i,j] ~ normal ( 0, sqrt(tau^2 * theta[i, j]^2) );
      if (n[i, j] == 0) {

        target += log_sum_exp(bernoulli_lpmf(1 | omega[j]),
        bernoulli_lpmf(0 | omega[j])
        + poisson_log_lpmf(0 | log_mu[i, j] ) );

      } else {

        target += bernoulli_lpmf(0 | omega[j])
        + poisson_log_lpmf(n[i, j] | log_mu[i, j] );

      }
    }

    beta_AE[i] ~ normal (0, sigma_AE);

  }

}

generated quantities {

  array[I, J] real<lower=0, upper=1> zi;
  array[I, J] int<lower=0> n_pred;
  array[I, J] int<lower=0, upper=1> zi_pred;
  array[I, J] real<lower=0> lambda;
  array[I, J] real<lower=0> lambda_tilde;

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
      lambda[i, j] = (1 - zi_pred[i, j]) * exp(log_lambda_tilde_total[i, j]);
      lambda_tilde[i, j] = exp(log_lambda_tilde_total[i, j]);
    }
  }

}

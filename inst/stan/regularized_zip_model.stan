data {

  int<lower=0> I;
  int<lower=0> J;
  array[I, J] real<lower=0> E;
  array[I, J] int<lower=0> n;

  real <lower =0> scale_global ;
  real <lower =1> nu_global; // degrees of freedom for the half -t prior for tau
  real <lower =1> nu_local; // degrees of freedom for the half -t priors for thetas
  real <lower =0> slab_scale; // slab scale for the regularized horseshoe
  real <lower =0> slab_df; // slab degrees of freedom for the regularized horseshoe

}

transformed data {

  array[I, J] real log_E;
  log_E = log(E);

}

parameters {

  array[I, J] real<lower=0> lambda;
  real<lower=0> tau; // global shrinkage parameter
  array[I, J] real<lower=0> theta; // local shrinkage parameter
  array[J] real<lower=0, upper=1> omega;
  real<lower=0> caux;

}

transformed parameters {

  array[I, J] real log_lambda;
  real<lower=0> c;

  log_lambda = log(lambda);
  c = slab_scale * sqrt(caux);

}

model {

  array[I, J] real theta_tilde;

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

  //tau ~ student_t(nu_global, 0, scale_global);
  tau ~ uniform(0.5, 2);
  omega ~ beta(0.5, 0.5); //Jeffrey's prior
  caux ~ inv_gamma (0.5* slab_df , 0.5* slab_df );

  for (j in 1 : J){
    for (i in 1 : I){
      theta[i,j] ~ student_t(nu_local , 0, 1);
      //theta_tilde[i, j] = sqrt(c^2 * theta[i, j]^2 / (c^2 + tau^2 * theta[i, j]));
      // log_lambda[i, j] ~ normal ( 0, tau * theta_tilde[i, j] );
      //log_lambda[i, j] ~ normal ( 0, tau * theta[i, j] );
      //log_lambda[i, j] ~ normal ( 0, 1);
      log_lambda[i, j] ~ normal ( 0, tau);
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

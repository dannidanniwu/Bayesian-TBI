data {
  int<lower=0> N;           // Total number of observations
  int<lower=1> D;           // Total number of covariates
  row_vector[D] X[N];       // Covariates
  int<lower=0> J;           // Number of treatment
  int<lower=0> K_bspline;   // Number of B-spline basis functions
  matrix[N, K_bspline] B[D]; // B-spline basis matrix for each observation for each covariate
  int<lower=0,upper=1>  Y[N];     // Outcome variable
  int<lower=1, upper=J> A[N]; // Treatment assignment for each observation
  int<lower=0> N_new;           // Number of new observations
  row_vector[D] Xnew[N_new];    // New covariates
  matrix[N_new, K_bspline] Bnew[D]; // New B-spline basis matrix for each observation for each covariate
  int<lower=1, upper=J> Anew[N_new]; // New treatment assignments for each observation

}

parameters {
  matrix[J, K_bspline] beta_spline[D];   // Coefficients for B-spline basis for each cluster for each covariate
  vector[D] m;
}

model {
  real yhat[N];
  for (d in 1:D) {
    for (j in 1:J) {
      beta_spline[d][j] ~ normal(0, 5);
    }
  }
  
  m ~ normal(0, 5);

  // Likelihood
  for (i in 1:N) {
    yhat[i] = dot_product(X[i], m);
    for (d in 1:D) {
      yhat[i] += dot_product(B[d][i], beta_spline[d][A[i]]);
    }
    Y[i] ~ bernoulli_logit(yhat[i]);
  }
}

generated quantities {
  real eta[N_new];
  vector[N] contrast_distr;
  
  for (i in 1:N_new) {
    eta[i] = dot_product(Xnew[i], m);
    contrast_distr[i] = 0;
    for (d in 1:D) {
      eta[i] += dot_product(Bnew[d][i], beta_spline[d][Anew[i]]);
      contrast_distr[i] += dot_product(B[d][i], beta_spline[d][2]) - dot_product(B[d][i], beta_spline[d][1]);
    }
  }
  
}
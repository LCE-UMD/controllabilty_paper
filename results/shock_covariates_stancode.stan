// generated with brms 2.12.0
functions {

  /* compute the logm1 link 
   * Args: 
   *   p: a positive scalar
   * Returns: 
   *   a scalar in (-Inf, Inf)
   */ 
   real logm1(real y) { 
     return log(y - 1);
   }
  /* compute the inverse of the logm1 link 
   * Args: 
   *   y: a scalar in (-Inf, Inf)
   * Returns: 
   *   a positive scalar
   */ 
   real expp1(real y) { 
     return exp(y) + 1;
   }
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  vector[N] Z_2_2;
  vector[N] Z_2_3;
  vector[N] Z_2_4;
  vector[N] Z_2_5;
  vector[N] Z_2_6;
  int<lower=1> NC_2;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // residual SD
  real<lower=1> nu;  // degrees of freedom or shape
  // parameters for student-t distributed group-level effects
  real<lower=1> df_1;
  vector<lower=0>[N_1] udf_1;
  real<lower=1> df_2;
  vector<lower=0>[N_2] udf_2;
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  matrix[M_2, N_2] z_2;  // standardized group-level effects
  cholesky_factor_corr[M_2] L_2;  // cholesky factor of correlation matrix
}
transformed parameters {
  vector[N_1] dfm_1;
  vector[N_2] dfm_2;
  vector[N_1] r_1_1;  // actual group-level effects
  matrix[N_2, M_2] r_2;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_2] r_2_1;
  vector[N_2] r_2_2;
  vector[N_2] r_2_3;
  vector[N_2] r_2_4;
  vector[N_2] r_2_5;
  vector[N_2] r_2_6;
  dfm_1 = sqrt(df_1 * udf_1);
  dfm_2 = sqrt(df_2 * udf_2);
  r_1_1 = dfm_1 .* (sd_1[1] * (z_1[1]));
  // compute actual group-level effects
  r_2 = rep_matrix(dfm_2, M_2) .* (diag_pre_multiply(sd_2, L_2) * z_2)';
  r_2_1 = r_2[, 1];
  r_2_2 = r_2[, 2];
  r_2_3 = r_2[, 3];
  r_2_4 = r_2[, 4];
  r_2_5 = r_2[, 5];
  r_2_6 = r_2[, 6];
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + Xc * b;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_2_2[J_2[n]] * Z_2_2[n] + r_2_3[J_2[n]] * Z_2_3[n] + r_2_4[J_2[n]] * Z_2_4[n] + r_2_5[J_2[n]] * Z_2_5[n] + r_2_6[J_2[n]] * Z_2_6[n];
  }
  // priors including all constants
  target += student_t_lpdf(b | 3,0,10);
  target += student_t_lpdf(Intercept | 3, 0, 10);
  target += student_t_lpdf(sigma | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += gamma_lpdf(nu | 3.325,0.1)
    - 1 * gamma_lccdf(1 | 3.325,0.1);
  target += gamma_lpdf(df_1 | 3.325,0.1);
  target += inv_chi_square_lpdf(udf_1 | df_1);
  target += gamma_lpdf(df_2 | 3.325,0.1);
  target += inv_chi_square_lpdf(udf_2 | df_2);
  target += student_t_lpdf(sd_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_1[1] | 0, 1);
  target += student_t_lpdf(sd_2 | 3, 0, 10)
    - 6 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(to_vector(z_2) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_2 | 2);
  // likelihood including all constants
  if (!prior_only) {
    target += student_t_lpdf(Y | nu, mu, sigma);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // compute group-level correlations
  corr_matrix[M_2] Cor_2 = multiply_lower_tri_self_transpose(L_2);
  vector<lower=-1,upper=1>[NC_2] cor_2;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_2) {
    for (j in 1:(k - 1)) {
      cor_2[choose(k - 1, 2) + j] = Cor_2[j, k];
    }
  }
}


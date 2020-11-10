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
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
}
transformed parameters {
  vector[N_1] dfm_1;
  vector[N_1] r_1_1;  // actual group-level effects
  dfm_1 = sqrt(df_1 * udf_1);
  r_1_1 = dfm_1 .* (sd_1[1] * (z_1[1]));
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + Xc * b;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
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
  target += student_t_lpdf(sd_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_1[1] | 0, 1);
  // likelihood including all constants
  if (!prior_only) {
    target += student_t_lpdf(Y | nu, mu, sigma);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}


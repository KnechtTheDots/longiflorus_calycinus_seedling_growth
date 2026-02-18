data{
  int n_long; // total number of observations (individuals x time)
  int n_short; // number of individuals
  int n_pred;
  vector[n_long] area; // seedling area at each observation
  vector[n_long] age; // seedling age at each observation
  array[n_long] int id_long; // seedling id
  array[n_short] int survive; // did the seedling survive until the 4th day?
  array[n_short] int id_short; // id of seedlings in the short data
  vector[n_short] seed; // seed size of each individual
  vector[n_short] final_size; // size of seedling before drought treatment
  vector[n_short] age_max; // age at the last time point, for the post-predictive
  vector[n_pred] size_pred;
  vector[n_pred] rgr_pred;
}
transformed data{
  // log of area (exponential growth)
  vector[n_long] log_area = log(area); 
  // standardized seed size for initial size submodel
  vector[n_short] seed_std = (seed - mean(seed))/sd(seed); 
  real mu_final_size = mean(final_size);
  real sd_final_size = sd(final_size);
  // standardize final size for the logistic regression
  vector[n_short] final_size_std = (final_size - mu_final_size)/sd_final_size;
}
parameters{
  matrix[2, n_short] Z; // relative growth rate of each individual
  real gamma; // intercept of the seed size -> log_size_0 sub model
  real beta; // slope of the seed size -> log_size_0 sub model
  vector<lower=0>[2] tau; // scales of rgr and log_size_0
  real<lower=0> sig_m; // scale of log size
  real<lower=0> rgr_bar; // mean relative growth rate
  cholesky_factor_corr[2] L_Omega; // cholesky factor of the correlation matrix between rgr and log_size_0
  real alpha_survive; // intercept of the survival logistic regression
  real beta_survive; // slope of the survival logistic regression
  real b_rgr_surv;
}
transformed parameters{
  vector[n_short] log_size_0; // log(size_0) 
  vector[n_short] rgr; // relative growth rate
  {
  matrix[2, n_short] Mu; // define matrix of means
  matrix[2, n_short] Beta; // matrix carrying the log_size_0 and relative growth rate
  for(i in 1:n_short){
    Mu[1, i] = rgr_bar; 
    Mu[2, i] = gamma + beta*seed_std[i];
  }
  Beta = Mu + diag_pre_multiply(tau, L_Omega) * Z;
  rgr = Beta[1,]';
  log_size_0 = Beta[2,]';
  }
  
}
model{
  to_vector(Z) ~ normal(0, 1); // 95% ~ -2 <-> 2
  tau ~ exponential(1); // ~ .025 <-> 3.7
  L_Omega ~ lkj_corr_cholesky(2); // mildly regularizing toward no correlation
  beta ~ normal(0, 1); // ~ -2 <-> 2
  gamma ~ normal(0, 1); // ~ -2 <-> 2
  rgr_bar ~ normal(0, .25); // ~ -.5 <-> .5
  sig_m ~ exponential(1); // ~ .025 <-> 3.7
  alpha_survive ~ normal(0, 1);
  beta_survive ~ normal(0, 1);
  b_rgr_surv ~ normal(0, 1);
  
  // mu = E(log_area | log_size, rgr, age)
  vector[n_long] mu = log_size_0[id_long] + rgr[id_long] .* age;
  
  // likelihood for log area
  log_area ~ normal(mu, sig_m);
  
  // likelihood for survival
  survive ~ bernoulli_logit(alpha_survive + beta_survive * final_size_std + b_rgr_surv .* Z[1,]');
}
generated quantities{
  // recover the correlation matrix of log_size_0 and rgr from the cholesky decomposition
  // don't need it for now but keep it so I remember matrix[2,2] R = multiply_lower_tri_self_transpose(L_Omega);
  
  // posterior predictive of size at the final age
  array[n_short] real y_rep = normal_rng(log_size_0[id_short] + rgr[id_short] .* age_max, sig_m);
  
  // posterior of survival probability over log_size
  vector[n_pred] p_size_pred = inv_logit(alpha_survive + beta_survive * size_pred);
  // posterior of survival over relative growth rate at 3 ages and 3 seed sizes
  vector[n_pred] p_rgr_age11_pred = inv_logit(alpha_survive + beta_survive * (exp(gamma + rgr_pred*11) - mu_final_size)/sd_final_size);
  vector[n_pred] p_rgr_age12_pred = inv_logit(alpha_survive + beta_survive * (exp(gamma + rgr_pred*12) - mu_final_size)/sd_final_size);
  vector[n_pred] p_rgr_age13_pred = inv_logit(alpha_survive + beta_survive * (exp(gamma + rgr_pred*13) - mu_final_size)/sd_final_size);
  vector[n_pred] p_rgr_age14_pred = inv_logit(alpha_survive + beta_survive * (exp(gamma + rgr_pred*14) - mu_final_size)/sd_final_size);
  vector[n_pred] p_rgr_minus_1sd_seed = inv_logit(alpha_survive + beta_survive * (exp(gamma - beta + rgr_pred*14) - mu_final_size)/sd_final_size);
  vector[n_pred] p_rgr_plus_1sd_seed = inv_logit(alpha_survive + beta_survive * (exp(gamma + beta + rgr_pred*14) - mu_final_size)/sd_final_size);
  
}

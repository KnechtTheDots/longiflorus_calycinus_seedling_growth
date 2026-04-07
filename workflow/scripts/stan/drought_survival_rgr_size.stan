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
  vector[n_pred] z_tilde; // vector of standardized x values to predict over
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
  vector[n_short] age_max_std = (age_max - mean(age_max))/sd(age_max);
}
parameters{
  vector[n_short] z_rgr; // relative growth rate of each individual
  vector[n_short] z_size_0;
  real gamma; // intercept of the seed size -> log_size_0 sub model
  real beta; // slope of the seed size -> log_size_0 sub model
  real<lower=0> tau; // scales of rgr
  real<lower=0> tau_size_0;
  real<lower=0> sig_m; // scale of log size
  real<lower=0> rgr_bar; // mean relative growth rate
  real alpha_survive; // intercept of the survival logistic regression
  real beta_rgr_survive; // slope of the survival logistic regression or rgr
  real beta_size_survive;
  cholesky_factor_corr[4] L_R; // cholesky factor for correlations between rgr, final age and seed size
 }
transformed parameters{
  vector[n_short] log_size_0 = gamma + beta*seed_std + z_size_0*tau_size_0; // log(size_0) 
  vector[n_short] rgr = rgr_bar + z_rgr*tau; // relative growth rate

 // likelihood for survival
  vector[n_short] logit_p = alpha_survive + beta_rgr_survive * z_rgr +
                                            beta_size_survive * final_size_std;
}
model{
  tau ~ exponential(1); // ~ .025 <-> 3.7
  tau_size_0 ~ exponential(1);
  beta ~ normal(0, 1); // ~ -2 <-> 2
  gamma ~ normal(0, 1); // ~ -2 <-> 2
  rgr_bar ~ normal(0, .25); // ~ -.5 <-> .5
  sig_m ~ exponential(1); // ~ .025 <-> 3.7
  alpha_survive ~ normal(0, 1);
  beta_rgr_survive ~ normal(0, 1);
  beta_size_survive ~ normal(0, 1);
  L_R ~ lkj_corr_cholesky(2);
  
  array[n_short] vector[4] phenos;
  for(i in 1:n_short){
    phenos[i] = [z_rgr[id_short[i]], age_max_std[id_short[i]], seed_std[id_short[i]], z_size_0[id_short[i]]]';
  }
  
  // mean = 0 and sd of each phenotype is on (they are standardized) so no
  // need to muliply the cholesky corr by the standard deviations because it 
  // would just return L_R anyway
  phenos ~ multi_normal_cholesky([0,0,0,0]',L_R);
  
  // mu = E(log_area | log_size, rgr, age)
  vector[n_long] mu = log_size_0[id_long] + rgr[id_long] .* age;
  
  // likelihood for log area
  log_area ~ normal(mu, sig_m);
  
  
  survive ~ bernoulli_logit(logit_p);
}
generated quantities{
  
  
  
  vector[n_short] log_lik;
  for(i in 1:n_short){
    log_lik[i] = bernoulli_logit_lpmf(survive[i] | alpha_survive + 
                                        beta_rgr_survive * z_rgr[i] +
                                        beta_size_survive * final_size_std[i]);
  }
  // recover the correlation matrix of log_size_0 and rgr from the cholesky decomposition
  // don't need it for now but keep it commented so I remember later if I want it
  //matrix[2,2] R = multiply_lower_tri_self_transpose(L_Omega);
  
  // posterior predictive of size at the final age
  array[n_short] real y_rep = normal_rng(log_size_0[id_short] + rgr[id_short] .* age_max, sig_m);
  vector[n_short] final_rep = exp(to_vector(y_rep));
  
  // posterior predictive of survival
  array[n_short] int surv_rep;
  {
    vector[n_short] z_size_rep = (exp(to_vector(y_rep)) - mu_final_size)/sd_final_size;
    surv_rep = bernoulli_logit_rng(alpha_survive + beta_rgr_survive * z_rgr[id_short] + 
                                                   beta_size_survive * z_size_rep[id_short]);
  }
  
  vector[n_short] w = to_vector(surv_rep)/mean(to_vector(surv_rep));
  real sigma_w = sd(w);
  
  vector[n_pred] p_pred;
  {
    // calculate the final size (standardized to the data's values) using the mean
    // standardized seed size (i.e. 0) and an age of 14
    vector[n_pred] rgr_pred = rgr_bar + tau*z_tilde;
    vector[n_pred] size_pred = exp(gamma + rgr_pred*14);
    vector[n_pred] z_size_pred = (size_pred - mu_final_size)/sd_final_size;
    p_pred = inv_logit(alpha_survive + beta_rgr_survive * z_tilde + 
                                       beta_size_survive * z_size_pred);
  }
  
  
  
  
  vector[n_pred] rgr_pred = rgr_bar + z_tilde * tau;
  
  
  matrix[4,4] R_phenos = multiply_lower_tri_self_transpose(L_R);
}

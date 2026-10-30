// scale_ssm_model.stan

data {
  int<lower=1> C;                                  
  int<lower=1> T;                                  
  int<lower=1> N_obs;                              
  array[N_obs] int<lower=1, upper=C> core_idx;     
  array[N_obs] int<lower=1, upper=T> grid_idx;     
  vector[N_obs] y_log;                            
  array[N_obs] int<lower=0, upper=1> obs_indic;    
  vector[C] L_det;                                 
  array[C] int<lower=1, upper=T> t_lo;             
  array[C] int<lower=1, upper=T> t_hi;             
  real<lower=0> lkj_eta;                           
  int<lower=1> K;                                  
  array[K] int<lower=1, upper=C> pair_i;           
  array[K] int<lower=1, upper=C> pair_j;          
}

parameters { 
// core-specific parameters
  vector[C] mu;                                    
  vector<lower=-0.99, upper=0.99>[C] rho;          
  vector[C] log_sigma_proc_raw;                    
  vector[C] log_tau_obs_raw;                       
// correlation structure
  cholesky_factor_corr[C] L_R;                    
  real z_mean;                                   
  real<lower=0> z_sd;                            
  vector[C] z_init_state;                          
  matrix[T, C] eta;                                
}

transformed parameters {
  matrix[C, C] R = 
  multiply_lower_tri_self_transpose(L_R); 
vector<lower=0>[C] sigma_proc = exp(log_sigma_proc_raw);
vector<lower=0>[C] tau_obs    = exp(log_tau_obs_raw);
  matrix[C, C] L_Sigma;
  for (i in 1:C) {
    for (j in 1:C) {
      L_Sigma[i, j] = sigma_proc[i] * L_R[i, j];
    }
  }
  matrix[T, C] eps;
  for (t in 1:T) {
    eps[t] = (L_Sigma * eta[t]')';
  }
  vector[C] x_init;
for (c in 1:C) {
    x_init[c] = mu[c] + sigma_proc[c] * z_init_state[c];
}
  matrix[T, C] x;
  for (t in 1:T)
    for (c in 1:C)
      x[t, c] = mu[c];
  for (c in 1:C) {
    x[t_lo[c], c] = x_init[c];
    for (t in (t_lo[c] + 1):t_hi[c]) {
      x[t, c] = mu[c] + rho[c] * (x[t - 1, c] - mu[c]) + eps[t, c];
    }
  }
  vector[K] z_offdiag;
  for (k in 1:K) {
   real r_k = fmax(-0.999, fmin(0.999, R[pair_i[k], pair_j[k]]));
   z_offdiag[k] = 0.5 * log((1 + r_k) / (1 - r_k));   
 }
}

model {
// Priors
  mu          ~ normal(2, 1);
  rho ~ normal(0.5, 0.3);
  log_sigma_proc_raw  ~ normal(log(0.4), .3);   
  log_tau_obs_raw     ~ normal(log(0.5), .5); 
  L_R         ~ lkj_corr_cholesky(lkj_eta);        
  z_init_state ~ std_normal();
  to_vector(eta) ~ std_normal();
// Hierarchical layer on the cross-core correlations
  z_mean ~ normal(0, 1);
  z_sd   ~ exponential(2);                      
  target += normal_lpdf(z_offdiag | z_mean, z_sd);
// Observation model 
  for (i in 1:N_obs) {
    real mean_obs = x[grid_idx[i], core_idx[i]];
    real sd_obs   = tau_obs[core_idx[i]];
    if (obs_indic[i] == 1) {
      y_log[i] ~ normal(mean_obs, sd_obs);
    } else {
// Detection limit
      target += normal_lcdf(L_det[core_idx[i]] | mean_obs, sd_obs);
    }
  }
}

generated quantities {
  vector[K] R_offdiag;
  for (k in 1:K) {
    R_offdiag[k] = R[pair_i[k], pair_j[k]];
  }
 real median_r = tanh(z_mean);
}

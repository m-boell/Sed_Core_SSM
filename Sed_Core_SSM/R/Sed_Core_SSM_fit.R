# Sed_Core_SSM_fit.R ####

# Header ####
# Requires cmdstanr and a working cmdstan installation. To install:
#   install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev",getOption("repos")))
#   cmdstanr::install_cmdstan(cores = 4)
# Falls back to rstan if cmdstanr is unavailable.

# Libraries ####
suppressPackageStartupMessages({
  library(readxl) 
  library(dplyr) 
  library(tidyr) 
  library(tibble) 
  library(posterior) 
  library(here)
})

# Path settings ####
STAN_FILE <- here("Stan", "Sed_Core_SSM_model.stan")
OUT_DIR   <- here("Stan_run")
dir.create(OUT_DIR, showWarnings = FALSE)

# Time grid settings ####
GRID_STEP <- 10
GRID_START <- 385
GRID_END <- 1795

# MCMC settings ####
LKJ_ETA <- 1.0 
SEED <- 20260513
NCHAINS <- 4
N_WARMUP <- 2000
N_SAMPLE <- 2000
ADAPT_DELTA <- 0.999  
MAX_TREEDEPTH <- 12 

set.seed(SEED)

# Data ####
file_info <- tibble::tribble(
  ~core_ID, ~file, ~skip, 
  "SMB_KC","Santa_Monica_Basin.xlsx", 4L,
  "SBB_SH2015", "SH2015_Santa_Barbara_Basin.xlsx", 5L,
  "SBB_B1992_214", "Baumgartner1992_Santa_Barbara_Basin_214.xlsx",5L,
  "SBB_B1992_224", "Baumgartner1992_Santa_Barbara_Basin_224.xlsx",5L
)

read_core <- function(file, skip) {
  d <- read_excel(here("Data", file), skip = skip,
                  col_names = TRUE, na = c("NA", "na", ""))
  d %>%
    select(core_ID, bin_median_year_AD, bin_start_year_AD,bin_end_year_AD, raw_SDR) %>%
    mutate(across(c(bin_median_year_AD, bin_start_year_AD, bin_end_year_AD, raw_SDR),
                  ~ suppressWarnings(as.numeric(.x)))) %>%
    filter(!is.na(core_ID), !is.na(bin_median_year_AD))
}

dat <- bind_rows(lapply(seq_len(nrow(file_info)), function(i) {
  with(file_info[i, ], read_core(file, skip))
}))

# Build time grid
core_levels <- file_info$core_ID
C <- length(core_levels)
grid_centers <- seq(GRID_START, GRID_END, by = GRID_STEP)
T_grid <- length(grid_centers)

# Computes the time bin into which an observation falls
dat <- dat %>%
  mutate(grid_idx = as.integer(
    round((bin_median_year_AD - GRID_START) / GRID_STEP)) + 1L, 
    grid_idx = pmin.int(pmax.int(grid_idx, 1L), T_grid),
    core_idx = match(core_ID, core_levels)) 
dat_used <- dat %>% filter(!is.na(raw_SDR))
N_obs    <- nrow(dat_used)

# Cores' temporal range
t_c_lo <- vapply(seq_len(C),
                 function(c) as.integer(min(dat_used$grid_idx[dat_used$core_idx == c])),
                 integer(1))
t_c_hi <- vapply(seq_len(C),
                 function(c) as.integer(max(dat_used$grid_idx[dat_used$core_idx == c])),
                 integer(1))

# Per-core detection limit = log(minimum non-zero SDR/2)
L_per_core <- vapply(seq_len(C), function(c) {
  v <- dat_used$raw_SDR[dat_used$core_idx == c & dat_used$raw_SDR > 0]
  if (length(v) == 0) -10 else log(min(v) / 2)
}, numeric(1))

dat_used <- dat_used %>%
  mutate(obs_indic = as.integer(raw_SDR > 0),
         y_log     = ifelse(raw_SDR > 0, log(raw_SDR), 0.0))   

# Core pairs
pair_grid <- expand.grid(i = seq_len(C), j = seq_len(C)) %>%
  filter(i < j) %>% arrange(i, j)
K       <- nrow(pair_grid) 
pair_i  <- pair_grid$i
pair_j  <- pair_grid$j

cat(sprintf("Data: %d obs across %d cores; %d censored zeros.\n",
            N_obs, C, sum(dat_used$obs_indic == 0)))
cat("Per-core grid ranges (years AD):\n")
for (c in seq_len(C))
  cat(sprintf("  %-18s  [%4d, %4d]\n",
              core_levels[c], grid_centers[t_c_lo[c]], grid_centers[t_c_hi[c]]))
cat("Detection limits (log scale):\n"); print(setNames(round(L_per_core, 3), core_levels))

# Stan data ####
stan_data <- list(
  C         = C,
  T         = T_grid,
  N_obs     = N_obs,
  core_idx  = dat_used$core_idx,
  grid_idx  = dat_used$grid_idx,
  y_log     = dat_used$y_log,
  obs_indic = dat_used$obs_indic,
  L_det     = L_per_core,
  t_lo      = t_c_lo,
  t_hi      = t_c_hi,
  lkj_eta   = LKJ_ETA,
  K         = K,
  pair_i    = pair_i,
  pair_j    = pair_j
)

saveRDS(stan_data,    file.path(OUT_DIR, "stan_data.rds"))
saveRDS(dat_used,     file.path(OUT_DIR, "dat_used.rds"))
saveRDS(grid_centers, file.path(OUT_DIR, "grid_centers.rds"))
saveRDS(core_levels,  file.path(OUT_DIR, "core_levels.rds"))
saveRDS(pair_grid,    file.path(OUT_DIR, "pair_grid.rds"))

# Fit model ####
use_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE)

if (use_cmdstan) {
  library(cmdstanr)
  mod <- cmdstan_model(STAN_FILE) 
  
  fit <- mod$sample(
    data            = stan_data,
    chains          = NCHAINS,
    parallel_chains = NCHAINS, 
    iter_warmup     = N_WARMUP,
    iter_sampling   = N_SAMPLE,
    seed            = SEED,
    adapt_delta     = ADAPT_DELTA,
    max_treedepth   = MAX_TREEDEPTH,
    refresh         = 200 
  )
  
# Save draws
fit$save_object(file = file.path(OUT_DIR, "fit_cmdstan.rds"))
draws <- fit$draws(format = "draws_array")
saveRDS(draws, file.path(OUT_DIR, "draws.rds"))

# Diagnostics and summaries ####
print(fit$diagnostic_summary())
print(fit$summary(variables = c("mu", "rho", "sigma_proc", "tau_obs", "z_mean", "z_sd", "median_r", "R_offdiag")), n = 50)

draws_array <- if (use_cmdstan) draws else
  posterior::as_draws_array(rstan::extract(fit, permuted = FALSE))

summary_tbl <- posterior::summarise_draws(
  posterior::subset_draws(draws_array,
                          variable = c("mu", "rho", "sigma_proc", "tau_obs", "z_mean", "z_sd","median_r", "R_offdiag")),
  mean, median, sd, ~quantile(.x, c(0.05, 0.95)), default_convergence_measures()
)
write.csv(summary_tbl, file.path(OUT_DIR, "Posterior_summary.csv"), row.names = FALSE)

# Rstan alternative ####
} else if (requireNamespace("rstan", quietly = TRUE)) {
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = NCHAINS)
  
  fit <- stan(file    = STAN_FILE,
              data    = stan_data,
              chains  = NCHAINS,
              warmup  = N_WARMUP,
              iter    = N_WARMUP + N_SAMPLE,
              seed    = SEED,
              control = list(adapt_delta   = ADAPT_DELTA,
                             max_treedepth = MAX_TREEDEPTH))
  
  saveRDS(fit, file.path(OUT_DIR, "fit_rstan.rds"))
  cat("\n--- Summary (key parameters) ---\n")
  print(fit, pars = c("mu", "rho", "sigma_proc", "tau_obs", "z_mean", "z_sd", "median_r","R_offdiag"),
        probs = c(0.05, 0.5, 0.95))
} else {
  stop("Neither cmdstanr nor rstan is installed. See header for installation hints.")
}

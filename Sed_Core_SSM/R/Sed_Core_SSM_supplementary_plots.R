# Sed_Core_SSM_supplementary_plots.R ####

# Libraries ####
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(posterior)
  library(readxl)
  library(here)
})

# Path settings ####
OUT_DIR <- here("Stan_run")
DATA_DIR <- here("Data")
FIG_DIR <- file.path(OUT_DIR, "Figures")
SUPP_DIR <- file.path(FIG_DIR, "Supplementary figures")
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

# Load fit ####
core_levels <- readRDS(file.path(OUT_DIR, "core_levels.rds"))
grid_centers <- readRDS(file.path(OUT_DIR, "grid_centers.rds"))
pair_grid <- readRDS(file.path(OUT_DIR, "pair_grid.rds"))
stan_data <- readRDS(file.path(OUT_DIR, "stan_data.rds"))

C <- length(core_levels)
T_grid <- length(grid_centers)
t_c_lo <- stan_data$t_lo
t_c_hi <- stan_data$t_hi

cmdstan_path <- file.path(OUT_DIR, "fit_cmdstan.rds")
rstan_path <- file.path(OUT_DIR, "fit_rstan.rds")

if (file.exists(cmdstan_path)) {
  fit <- readRDS(cmdstan_path)
  draws_array <- fit$draws(format = "draws_array")
} else if (file.exists(rstan_path)) {
  fit <- readRDS(rstan_path)
  draws_array <- posterior::as_draws_array(rstan::extract(fit, permuted = FALSE))
} else {
  stop("No fit object found in stan_run/. Run SCB_Sed_Core_SSM_fit.R first.")
}

post_array <- function(name) {
  posterior::subset_draws(draws_array, variable = name) |>
    posterior::merge_chains() |>
    posterior::as_draws_matrix()
}

# Plot theme ####
theme_report <- theme_classic(base_size = 12) + theme(
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "grey92", color = NA),
  strip.text = element_text(face = "bold", size = 11.3),
  plot.title = element_text(face = "bold", size = 13.3)
)

# Supplementary Figure 1 ####
# Pairwise correlation posteriors
R_off <- post_array("R_offdiag")
pair_lab <- with(pair_grid, paste(core_levels[i], core_levels[j], sep = " - "))

pair_df <- tibble(pair = factor(rep(pair_lab, each = nrow(R_off)), levels = pair_lab),
                  r = as.numeric(R_off))

# % of posterior on either side of 0
pair_pct <- pair_df |>
  group_by(pair) |>
  summarise(p_neg = mean(r < 0), p_pos = mean(r > 0), med = median(r)) |>
  mutate(lab = sprintf("< 0: %.0f%%\n>0: %.0f%%", 100*p_neg, 100*p_pos))
print(pair_pct, n = Inf)

p_pairs <- ggplot(pair_df, aes(r, fill = r > 0)) +
  geom_histogram(binwidth = 0.035, boundary = 0, color = NA) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_vline(data = pair_pct, aes(xintercept = med), color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_text(data = pair_pct, aes(label = lab), x = Inf, y = Inf,
            hjust = 1.1, vjust = 1.3, size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"), guide = "none") +
  facet_wrap(~ pair, ncol = 3, scales = "free_y") +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = "correlation", y = "posterior density",
       title = "Pairwise correlation posteriors") +
  theme_report
print(p_pairs)
ggsave(file.path(SUPP_DIR, "Supplementary_Fig_1.png"), p_pairs, width = 11, height = 7, dpi = 600)

# Supplementary Figure 2 ####
# SDR & covariate records
# Plot settings
XLIM       <- c(950, 1450)
MAJOR_STEP <- 50  
MINOR_STEP <- 10   
CRI_COLS   <- c("q25", "q75")  
SHOW_ALL_X_LABELS <- TRUE 
FONT_MIN <- 12
FONT_MAX <- 14
FIG_W <- 9; FIG_H <- 12; FIG_DPI <- 600

# A. SDR posteriors
x_post  <- post_array("x")
n_draws <- nrow(x_post)
x_arr <- array(NA_real_, dim = c(n_draws, T_grid, C))
for (ci in seq_len(C)) {
  for (ti in seq_len(T_grid)) {
    col_name <- sprintf("x[%d,%d]", ti, ci)
    if (col_name %in% colnames(x_post)) x_arr[, ti, ci] <- x_post[, col_name]
  }
}
x_summary <- expand.grid(t = seq_len(T_grid), c = seq_len(C)) |>
  mutate(year = grid_centers[t],
         core = core_levels[c],
         med  = mapply(function(t, c) median(x_arr[, t, c]), t, c),
         q05  = mapply(function(t, c) quantile(x_arr[, t, c], 0.05), t, c),
         q25  = mapply(function(t, c) quantile(x_arr[, t, c], 0.25), t, c),
         q75  = mapply(function(t, c) quantile(x_arr[, t, c], 0.75), t, c),
         q95  = mapply(function(t, c) quantile(x_arr[, t, c], 0.95), t, c)) |>
  mutate(core = factor(core, levels = core_levels))

# Core ranges
range_df <- tibble(core = core_levels,
                   year_lo = grid_centers[t_c_lo],
                   year_hi = grid_centers[t_c_hi]) |>
  mutate(core = factor(core, levels = core_levels))
x_summary <- x_summary |>
  left_join(range_df, by = "core") |>
  filter(year >= year_lo & year <= year_hi)

# Z-score the records
x_scaled <- x_summary |>
  group_by(core) |>
  mutate(core_mid    = mean(med, na.rm = TRUE),
         core_spread = sd(med, na.rm = TRUE),
         med_s = (med - core_mid) / core_spread,
         lo_s  = (.data[[CRI_COLS[1]]] - core_mid) / core_spread,
         hi_s  = (.data[[CRI_COLS[2]]] - core_mid) / core_spread) |>
  ungroup()

# B. SST
sst <- read_excel(file.path(DATA_DIR, "McCabe_Glynn_et_al_2013_Speleothem_O_isotopes.xlsx")) |>
  transmute(year = as.numeric(age_AD), sst = as.numeric(sst.anom)) |>
  filter(!is.na(year), !is.na(sst)) |>
  arrange(year) |>
  mutate(sst_z = as.numeric(scale(sst)))

# C. Mo
mo <- read_excel(file.path(DATA_DIR, "Wang_Hendy_2021_Elemental_Data_SBB.xlsx")) |>
  transmute(year = as.numeric(age_CE), Mo = as.numeric(Mo)) |>
  filter(!is.na(year), !is.na(Mo)) |>
  arrange(year) |>
  mutate(Mo = as.vector(scale(Mo)))

# D. Upwelling
rhiz <- read_excel(file.path(DATA_DIR, "Barron_et_al_2015_Diatom_Data.xlsx"), sheet = 1) |>
  filter(grepl("^Rhizosolenia", Species)) |>
  transmute(year = as.numeric(Year), rhiz = as.numeric(Relative_Percent_Abundance)) |>
  filter(!is.na(year), !is.na(rhiz)) |>
  arrange(year) |>
  mutate(rhiz = as.vector(scale(rhiz)))

# More plot settings
core_colors <- c(SMB_KC        = "#000000",
                 SBB_SH2015    = "#E69F00",
                 SBB_B1992_214 = "#56B4E9",
                 SBB_B1992_224 = "#CC79A7")
record_colors <- c(core_colors[core_levels],
                   sst  = "#999999",
                   mo   = "#009E73",
                   rhiz = "#0072B2")
SST_WARM <- "red"
SST_COOL <- "#2c7fb8"
record_labels <- c(setNames(toupper(core_levels), core_levels),
                   sst = "SST", mo = "Mo", rhiz = "Rhizosolenia")

theme_stack <- function() {
  theme_classic(base_size = FONT_MIN) +
    theme(text                    = element_text(size = FONT_MIN),
          axis.text               = element_text(size = FONT_MIN, colour = "black"),
          axis.title              = element_text(size = FONT_MAX),
          legend.text             = element_text(size = FONT_MIN),
          legend.title            = element_text(size = FONT_MIN),
          plot.tag                = element_text(size = FONT_MAX, face = "bold"),
          axis.ticks.length       = unit(6, "pt"),
          axis.minor.ticks.length = rel(0.5),
          plot.margin             = margin(t = 5.5, r = 28, b = 5.5, l = 2.5, unit = "pt"))
}
major_breaks <- seq(ceiling(XLIM[1] / MAJOR_STEP) * MAJOR_STEP,
                    floor(XLIM[2]   / MAJOR_STEP) * MAJOR_STEP, by = MAJOR_STEP)
minor_breaks <- setdiff(seq(ceiling(XLIM[1] / MINOR_STEP) * MINOR_STEP,
                            floor(XLIM[2]   / MINOR_STEP) * MINOR_STEP, by = MINOR_STEP),
                        major_breaks)
x_scale <- function() {
  scale_x_continuous(breaks = major_breaks, minor_breaks = minor_breaks,
                     guide = guide_axis(minor.ticks = TRUE), expand = c(0, 0))
}
hide_x_labels <- if (SHOW_ALL_X_LABELS) theme() else theme(axis.text.x = element_blank())

# Panels
# A. SDRs
p_cores <- ggplot(x_scaled, aes(x = year, colour = core, fill = core)) +
  geom_hline(yintercept = 0, colour = "black", linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lo_s, ymax = hi_s), alpha = 0.15, colour = NA) +
  geom_line(aes(y = med_s), linewidth = 0.9) +
  scale_colour_manual(values = core_colors, labels = record_labels[core_levels], name = NULL) +
  scale_fill_manual(values = core_colors, guide = "none") +
  x_scale() +
  coord_cartesian(xlim = XLIM) +
  labs(x = NULL, y = "latent ln(SDR)") +
  theme_stack() + hide_x_labels +
  theme(legend.position = "top")

# B. SSTs
p_sst <- ggplot(sst, aes(x = year)) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(sst_z, 0)), fill = SST_WARM) +
  geom_ribbon(aes(ymin = pmin(sst_z, 0), ymax = 0), fill = SST_COOL) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = .8) +
  x_scale() +
  coord_cartesian(xlim = XLIM) +
  labs(x = NULL, y = expression("SST anomaly (" * degree * "C)")) +
  theme_stack() + hide_x_labels

# C. Mo (oxygen)
p_mo <- ggplot(mo, aes(x = year, y = Mo)) +
  geom_hline(yintercept = 0, colour = "#009E73", linetype = "dashed", linewidth = 0.8) +
  geom_line(colour = record_colors[["mo"]], linewidth = 1) +
  x_scale() +
  coord_cartesian(xlim = XLIM) +
  labs(x = NULL, y = expression(atop(Mo[EF], "(z-score)"))) +
  theme_stack()

# D. Upwelling
p_rhiz <- ggplot(rhiz, aes(x = year, y = rhiz)) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = "dashed", linewidth = 0.8) +
  geom_line(colour = record_colors[["rhiz"]], linewidth = 1) +
  x_scale() +
  labs(x = "Year AD", y = expression(atop(italic("Rhizosolenia") ~ "spp.", "(z-score)"))) +
  coord_cartesian(xlim = XLIM) +
  theme_stack()

# Final figure
fig_covariates <- (p_cores / p_sst / p_mo / p_rhiz) +
  plot_layout(heights = c(2.75, 2.75, 2.75, 2.75)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = FONT_MAX, face = "bold", hjust = 1.5, vjust = -0.10, margin = margin(r = 5, b = 5))))

ggsave(file.path(SUPP_DIR, "Supplementary_Fig_2.png"), fig_covariates, width = FIG_W, height = FIG_H, dpi = FIG_DPI)

# Supplementary Figure 3 ####
# SMB vs. SBB CalCOFI temperature comparison
# Settings
infile <- here("Data", "SMB_SBB_CalCOFI_temperature_data.xlsx")
min_bins <- 5             
yr_range <- c(1984, 2021) 
max_gap  <- 8

# Load data
d <- read_excel(infile, col_types = c("text", "text", "text", "numeric", "numeric", "numeric", "numeric"))
d <- d[complete.cases(d$SBB, d$SMB), ]
d$date <- as.Date(sprintf("%s-%s-15", d$Year, d$Month))

# Mean temp/cruise/basin
n  <- aggregate(SMB ~ date, d, length)
tt <- aggregate(cbind(SBB, SMB) ~ date, d, mean)
tt <- tt[n$SMB >= min_bins, ]
yr <- as.numeric(format(tt$date, "%Y"))
tt <- tt[yr >= yr_range[1] & yr <= yr_range[2], ]
tt <- tt[order(tt$date), ]

tt$seg <- cumsum(c(FALSE, diff(as.numeric(tt$date)) / 30.44 > max_gap))

long <- rbind(data.frame(tt[c("date", "seg")], basin = "SMB", temp = tt$SMB),
              data.frame(tt[c("date", "seg")], basin = "SBB", temp = tt$SBB))
long$basin <- factor(long$basin, levels = c("SMB", "SBB"))

# Plot
p <- ggplot(long, aes(date, temp, colour = basin, group = interaction(basin, seg))) +
  geom_line(linewidth = 0.6) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y")+
  scale_colour_manual(values = c(SMB = "#b2182b", SBB = "#2166ac")) +
  labs(title = "CalCOFI temperature comparison", subtitle = "SMB: line 86.7, station 35.0\nSBB: line 81.8, station 46.9", x = "Year", y = "Mean temperature, 0-100 m (\u00B0C)", colour = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", legend.justification = "center") +
  theme(title = element_text(face = "bold"))
ggsave(file.path(SUPP_DIR, "Supplementary_Fig_3.png"), p, width = 10, height = 4.5, dpi = 600)

# Supplementary Figure 4 ####
# Scalar forest plot
make_post_summary <- function(name) {
  m <- post_array(name)
  tibble(core = core_levels,
         med  = apply(m, 2, median),
         lo   = apply(m, 2, quantile, 0.05),
         hi   = apply(m, 2, quantile, 0.95)) |>
    mutate(core = factor(core, levels = rev(core_levels)))
}

scalar_summaries <- bind_rows(
  make_post_summary("mu")         |> mutate(param = "mu"),
  make_post_summary("rho")        |> mutate(param = "rho"),
  make_post_summary("sigma_proc") |> mutate(param = "sigma_proc"),
  make_post_summary("tau_obs")    |> mutate(param = "tau_obs")
) |> mutate(param = factor(param, levels = c("mu", "rho", "sigma_proc", "tau_obs")))

p_scalar <- ggplot(scalar_summaries, aes(x = med, y = core)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.18,
                color = "#2c7fb8", orientation = "y") +
  geom_point(color = "#08519c", size = 2) +
  facet_wrap(~ param, scales = "free_x", nrow = 1) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  labs(x = "posterior median (90% CrI)", y = NULL,
       title = "Scalar parameter posteriors") +
  theme_report

ggsave(file.path(SUPP_DIR, "Supplementary_Fig_4.png"), p_scalar, width = 11.5, height = 3.6, dpi = 600)
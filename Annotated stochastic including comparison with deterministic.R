############################################
# VAPOR: STOCHASTIC SINGLE-GENERATION MODEL
# Maps deterministic R₀ onto stochastic outcomes
# No CFD scaling applied yet (for single-gen only)
# Author: David Fisman et al.
# Updated: June 2025
############################################

library(MASS)
library(ggplot2)

# ---------------------------
# Core Parameters
# ---------------------------
N <- 100                 # Total susceptible population
q1 <- 0.00695            # Close-contact transmission probability
p <- 0.5                 # Pulmonary ventilation rate (m³/hour)
V <- 500                 # Room volume (m³)
gen_time <- 96           # Generation time in hours (e.g., 4 days)
initial_infected <- 1    # Index case count
num_runs <- 1000         # Number of stochastic simulations per parameter combo

# ---------------------------
# Deterministic R₀ Calculator
# ---------------------------
calc_deterministic_r0 <- function(q, ACH, f) {
  P_a <- 1 - exp(-q * p * gen_time / (V * ACH))           # Wells-Riley (airborne) risk
  R0_total <- N * ((1 - f) * q1 + f * (q1 + P_a))          # Sum of contact and aerosol components
  return(R0_total)
}

# ---------------------------
# Stochastic Simulation (1 Generation)
# ---------------------------
run_sim <- function(q, ACH, f) {
  S <- N - initial_infected
  C <- initial_infected
  
  # Partition index case into Cb (aerosolizing) and Ca (contact only)
  Cb <- rbinom(1, C, f)
  Ca <- C - Cb
  
  # Probabilities from each route
  P_c <- 1 - (1 - q1)^(Ca + Cb)                              # Close contact
  P_a <- 1 - exp(-Cb * q * p * gen_time / (V * ACH))        # Airborne (Wells-Riley)
  
  # Combined risk (independent probabilities)
  P_total <- 1 - (1 - P_c) * (1 - P_a)
  
  # Sample number of new cases
  return(rbinom(1, S, min(1, P_total)))
}

# ---------------------------
# Parameter Grid + Simulation
# ---------------------------
ACH_vals <- c(0.5, 1, 2, 3)                # Ventilation rates (ACH)
f_vals <- c(0.05, 0.1, 0.15)               # Fraction of aerosolizing index cases
q_vals <- seq(5, 100, by = 5)              # Quanta generation rates

results <- data.frame()

for (ACH in ACH_vals) {
  for (f in f_vals) {
    for (q in q_vals) {
      det_R0 <- calc_deterministic_r0(q, ACH, f)                # Deterministic baseline
      
      gen1_cases <- replicate(num_runs, run_sim(q, ACH, f))     # Stochastic realizations
      stoch_R0 <- mean(gen1_cases)                              # Mean secondary cases
      
      results <- rbind(results, data.frame(
        q = q,
        ACH = as.factor(ACH),
        f = as.factor(f),
        R0_det = det_R0,
        R0_stoch = stoch_R0
      ))
    }
  }
}

# ---------------------------
# Plot: Stochastic vs. Deterministic R₀
# ---------------------------
p2a <- ggplot(results, aes(x = R0_det, y = R0_stoch, color = f, shape = ACH)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  #labs(
   # title = expression("Stochastic vs. Deterministic " * R[0]),
    #x = expression("Deterministic " * R[0] * " (per index case)"),
    #y = expression("Stochastic " * R[0] * " (mean 1-gen cases)"),
    #color = "Aerosolizing fraction (f)",
    #shape = "ACH"
  #) +
  labs(
    x = expression("Deterministic " * R[0] * " (per index case)"),
    y = expression("Stochastic " * R[0] * " (mean 1-gen cases)"),
    color = "Aerosolizing fraction (f)",
    shape = "ACH"
  ) +
  theme_minimal(base_size = 14)

print(p2a)
############################################
# VAPOR: FULL STOCHASTIC MODEL
# Includes R0 and k estimation from Gen 1
# with parameter sensitivity grid and visualization
# Author: David Fisman et al.
# Updated: June 2025
############################################

# ---------------------------
# Libraries
# ---------------------------
library(MASS)
library(ggplot2)

# ---------------------------
# Core Parameters
# ---------------------------
N <- 100
q1 <- 0.00695
q <- 20
p <- 0.5
V <- 500
ACH <- 0.5
dt <- 96
frac_aerosol <- 0.1
kappa <- 0.07
initial_infected <- 1
num_trials <- 1000

# ---------------------------
# Full Stochastic Gen 1 Model
# ---------------------------
run_full_stochastic_gen1 <- function(N = 100, q1 = 0.00695, q = 20, p = 0.5,
                                     V = 500, ACH = 0.5, dt = 96,
                                     frac_aerosol = 0.1, kappa = 0.07, initial_infected = 1) {
  S <- N - initial_infected
  C <- initial_infected
  
  Cb <- rbinom(1, C, frac_aerosol)
  Ca <- C - Cb
  
  P_c <- 1 - (1 - q1)^(kappa * (Ca + Cb))
  P_a <- 1 - exp(-kappa * Cb * q * p * dt / (V * ACH))
  P_total <- 1 - (1 - P_c) * (1 - P_a)
  P_total <- min(1, P_total)
  
  new_C <- rbinom(1, S, P_total)
  return(new_C)
}

# ---------------------------
# Estimate R0 and k
# ---------------------------
estimate_params <- function(data) {
  R0 <- mean(data)
  if (length(unique(data)) > 1 && sum(data > 0) >= 10) {
    fit <- tryCatch(
      suppressWarnings(fitdistr(data, "Negative Binomial")),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      k <- fit$estimate["size"]
      return(list(R0 = R0, k = k, fit = fit))
    }
  }
  return(list(R0 = R0, k = NA, fit = NULL))
}

# ---------------------------
# Run Simulation and Estimate
# ---------------------------
set.seed(123)
gen1_cases <- replicate(num_trials, run_full_stochastic_gen1())
gen1_stats <- estimate_params(gen1_cases)

cat("Full-Stochastic R0:", round(gen1_stats$R0, 2), "  k:", round(gen1_stats$k, 2), "\n")

# ---------------------------
# Plot NB Fits and Histogram
# ---------------------------

library(ggplot2)
library(patchwork)  # install.packages("patchwork") if needed

# Vectors: your simulated data
x_vals <- 0:max(gen1_cases, 30)
sim_df <- data.frame(gen1_cases = gen1_cases)

# Prepare NB fit curves
ref_fits <- data.frame(
  Pathogen = c("VAPOR Model", "SARS-1 Beijing", "SARS-1 Singapore", "1918 Flu", "MERS-CoV"),
  R0 = c(gen1_stats$R0, 2, 3, 2, 0.7),
  k = c(gen1_stats$k, 0.16, 0.11, 0.94, 0.25),
  color = c("blue", "red", "darkorange", "forestgreen", "purple"),
  lty = c("solid", "dashed", "dotdash", "twodash", "longdash")
)

fit_df <- do.call(rbind, lapply(1:nrow(ref_fits), function(i) {
  data.frame(
    SecondaryCases = x_vals,
    Density = dnbinom(x_vals, size = ref_fits$k[i], mu = ref_fits$R0[i]),
    Pathogen = ref_fits$Pathogen[i]
  )
}))

# Panel A: Histogram Only
pA <- ggplot(sim_df, aes(x = gen1_cases)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.6) +
  #labs(title = "A) Simulated Gen 1 Secondary Infections",
      # x = "Secondary Cases", y = "Density") +
  labs(x = "Secondary Cases", y = "Density") +
  theme_minimal(base_size = 14)

# Panel B: NB Fits
pB <- ggplot(fit_df, aes(x = SecondaryCases, y = Density, color = Pathogen, linetype = Pathogen)) +
  geom_line(size = 1.3) +
  #labs(title = "B) NB Model Fits from Simulation & Literature",
   #    x = "Secondary Cases", y = "Density") +
  labs(x = "Secondary Cases", y = "Density") +
  scale_color_manual(values = setNames(ref_fits$color, ref_fits$Pathogen)) +
  scale_linetype_manual(values = setNames(ref_fits$lty, ref_fits$Pathogen)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

#Combine plots (side by side) with patchwork
(pA | pB) + plot_layout(widths = c(1, 1.2))

# ---------------------------
# Sensitivity Grid
# ---------------------------
param_grid <- expand.grid(
  q = c(25, 50, 75),
  kappa = c(0.05, 0.1, 0.15),
  frac_aerosol = c(0.01, 0.05, 0.1),
  ACH = c(0.5, 1, 2)
)

estimate_R0_grid <- function(q, kappa, frac_aerosol, ACH, runs = 10000) {
  gen1s <- replicate(runs, run_full_stochastic_gen1(q = q, kappa = kappa,
                                                    frac_aerosol = frac_aerosol,
                                                    ACH = ACH))
  mean(gen1s, na.rm = TRUE)
}

# Compute R0 over grid
results <- do.call(rbind, lapply(1:nrow(param_grid), function(i) {
  row <- param_grid[i, ]
  R0_est <- estimate_R0_grid(q = row$q, kappa = row$kappa,
                             frac_aerosol = row$frac_aerosol,
                             ACH = row$ACH)
  data.frame(row, R0 = R0_est)
}))

#---------------------------------------------------------------------------
#CREATE 3-PANEL VERSION OF FIGURE 2-----------------------------------------
#---------------------------------------------------------------------------

library(patchwork)

# Adjust margins if needed for tight layout
p2a <- p2a + theme(plot.margin = margin(10, 10, 10, 10))
p2b <- pA + theme(plot.margin = margin(10, 10, 10, 10))
p2c <- pB + theme(plot.margin = margin(10, 10, 10, 10))

# Compose horizontal layout: (a) | (b) and (c) stacked vertically
fig2 <- (p2a / (p2b | p2c)) +
  plot_annotation(
    #title = "Figure 2 | Stochastic secondary transmission and offspring distributions",
    tag_levels = 'a',
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      plot.tag = element_text(size = 14, face = "bold")
    )
  )

print(fig2)

# Save figure
ggsave("Figure2_Composite.png", fig2, width = 12, height = 6.5, dpi = 300)
ggsave("Figure2_Composite.pdf", fig2, width = 12, height = 6.5)


# ---------------------------
# Visualize Grid Results
# ---------------------------
library(ggplot2)

p_grid <- ggplot(results, aes(x = ACH, y = R0, color = factor(q))) +
  geom_line(aes(group = interaction(q, kappa, frac_aerosol))) +
  facet_grid(
    rows = vars(frac_aerosol),
    cols = vars(kappa),
    labeller = labeller(
      frac_aerosol = function(x) paste("Aerosol Prob.:", x),
      kappa = function(x) paste("kappa:", x)
    )
  ) +
#  labs(
#    title = expression(R[0] ~ "Sensitivity Grid (Stochastic Gen 1)"),
#    x = "Air Changes per Hour (ACH)",
#    y = expression(R[0]),
#    color = "Quanta Rate (q)"
#  ) +
  labs(
    x = "Air Changes per Hour (ACH)",
    y = expression(R[0]),
    color = "Quanta Rate (q)"
  ) +
  theme_minimal()

# Display in console
print(p_grid)

# Save as high-quality PNG
ggsave("R0_Sensitivity_Grid_Stochastic.png", p_grid, width = 12, height = 8, dpi = 300)
# Save as PDF
ggsave("R0_Sensitivity_Grid_Stochastic.pdf", p_grid, width = 12, height = 8)












# -------------------------------------------------
# BLOCK 1: Stochastic Patch Model Simulation
# -------------------------------------------------
# Simulates stochastic spread in a multi-patch system with patch-specific ventilation (ACH)
# and a customizable exposure/contact matrix.
#
# Args:
#   kappa: Air mixing parameter (unitless; proxy for incomplete mixing)
#   f: Fraction of infectives who aerosolize (Cb)
#   V: Room volume (m3)
#   q: Quanta generation rate (quanta/hour per aerosolizer)
#   p: Pulmonary ventilation rate (m3/hour)
#   q1: Close-contact transmission probability per generation
#   dt: Duration of a generation (hours)
#   ACH_vals: Vector of ACH for each patch (length = number of patches)
#   pop_sizes: Vector of population sizes per patch (length = number of patches)
#   exposure_matrix: Matrix [patch x patch] of contact intensities (default: mostly diagonal)
#   max_time: Number of generations to simulate (integer)
#   seed_patch: Which patch to seed with initial case (default: 1)
#   debug: If TRUE, prints detailed transmission probabilities for each patch/generation
#
# Returns:
#   Matrix of incident cases by patch and generation (rows = generations, cols = patches)

run_stochastic_patch <- function(
    kappa = 0.15, f = 0.1, 
    V = 500, q = 50, p = 0.5, q1 = 0.00695, dt = 96, 
    ACH_vals = c(0.5, 2, 2), pop_sizes = c(100, 100, 100),
    exposure_matrix = matrix(c(1, 0.01, 0.01,
                               0.01, 1, 0.01,
                               0.01, 0.01, 1), nrow = 3, byrow = TRUE),
    max_time = 30, seed_patch = 1, debug = FALSE
) {
  num_patches <- length(pop_sizes)
  S <- pop_sizes                  # Susceptibles in each patch
  C <- rep(0, num_patches)        # Infectious in each patch
  R <- rep(0, num_patches)        # Removed (recovered/immune) in each patch
  
  # Seed initial infection in the specified patch
  C[seed_patch] <- 1
  S[seed_patch] <- S[seed_patch] - 1
  
  # Track incidence over time (rows: time, cols: patches)
  incidence <- matrix(0, nrow = max_time + 1, ncol = num_patches)
  incidence[1, ] <- C
  
  for (t in 1:max_time) {
    new_C <- rep(0, num_patches)
    for (i in 1:num_patches) {
      if (C[i] == 0 || S[i] == 0) next # Skip if no infectious or susceptibles
      
      # Partition infectious in patch i into aerosolizers (Cb) and others (Ca)
      Cb <- rbinom(1, C[i], f)
      Ca <- C[i] - Cb
      
      for (j in 1:num_patches) {
        if (S[j] == 0) next # No susceptibles in patch j
        
        # Effective exposure from patch i to j
        exposure_factor <- exposure_matrix[i, j] * kappa
        
        # Close-contact risk
        P_c <- 1 - (1 - q1)^(exposure_factor * (Ca + Cb))
        
        # Airborne/aerosol risk
        P_a <- 1 - exp(-exposure_factor * Cb * q * p * dt / (V * ACH_vals[j]))
        
        # Combined risk (independent risks)
        P_total <- 1 - (1 - P_c) * (1 - P_a)
        P_total <- min(1, P_total)  # Bound probability at 1
        
        # Print details if debugging
        if (debug) {
          cat(sprintf("t=%d | From Patch %d to Patch %d | P_c=%.4f | P_a=%.4f | P_total=%.4f\n",
                      t, i, j, P_c, P_a, P_total))
        }
        
        # Draw number of new infections in patch j caused by i
        new_infections <- rbinom(1, S[j], P_total)
        new_C[j] <- new_C[j] + new_infections
        S[j] <- S[j] - new_infections
      }
    }
    R <- R + C      # Move current infectious to removed
    C <- new_C      # Update infectious for next time step
    incidence[t + 1, ] <- C
  }
  return(incidence)
}

# -------------------------------------------------
# BLOCK 2: Multi-Run Simulation and Tidy Results
# -------------------------------------------------
# Runs the stochastic patch model repeatedly and aggregates final outbreak sizes
# across all patches and all simulation runs.
#
# Args:
#   num_runs: Number of independent simulations to run
#   ...: Additional parameters passed to run_stochastic_patch()
#
# Returns:
#   Data frame: each row = a run × patch, columns: patch, final_size, run_id

run_multiple_patch_simulations <- function(num_runs = 1000, ...) {
  sim_list <- vector("list", num_runs)
  for (i in seq_len(num_runs)) {
    sim_mat <- run_stochastic_patch(...)
    # Final outbreak size in each patch (sum all generations)
    sim_list[[i]] <- colSums(sim_mat)
  }
  sim_matrix <- do.call(rbind, sim_list)
  patch_names <- paste0("Patch_", seq_len(ncol(sim_matrix)))
  df <- as.data.frame(sim_matrix)
  colnames(df) <- patch_names
  df$run_id <- seq_len(num_runs)
  # Reshape to tidy/long format for plotting
  library(tidyr)
  df_long <- df %>%
    pivot_longer(cols = starts_with("Patch_"),
                 names_to = "patch",
                 values_to = "final_size")
  return(df_long)
}

# Example usage
set.seed(123)
sim_results <- run_multiple_patch_simulations(
  num_runs = 10000,
  kappa = 0.15,
  f = 0.1,
  ACH_vals = c(0.5, 6, 6),
  seed_patch = 1,
  debug = FALSE
)

# -------------------------------------------------
# BLOCK 3: Visualization & Summary Statistics
# -------------------------------------------------
library(ggplot2)
library(dplyr)

# Boxplot & jitter: total outbreak size per patch
ggplot(sim_results, aes(x = patch, y = final_size,  color=patch)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.3) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("Patch_1" = "skyblue",
                                "Patch_2" = "salmon",
                                "Patch_3" = "lightgreen")) +
  labs(
    title = "Final Outbreak Size per Patch (1000 runs)",
    x = "Patch",
    y = "Total Cases"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# Violin plot as alternative (shows shape/distribution)
ggplot(sim_results, aes(x = patch, y = final_size, fill = patch)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.08, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Patch_1" = "skyblue",
                               "Patch_2" = "salmon",
                               "Patch_3" = "lightgreen")) +
  labs(
    title = "Final Outbreak Size per Patch (Violin)",
    x = "Patch",
    y = "Total Cases"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Summary statistics: mean, median, 95% CI, extinction
patch_summary <- sim_results %>%
  group_by(patch) %>%
  summarize(
    mean_size = mean(final_size),
    median_size = median(final_size),
    ci_lower = quantile(final_size, 0.025),
    ci_upper = quantile(final_size, 0.975),
    extinction_rate = mean(final_size <= 1),
    .groups = "drop"
  )

print(patch_summary)

#-------------------------------------------------
# BLOCK 4: Simulate Different Ventilation Scenarios (ACH)
# With scenario summary columns and stratified plotting
#-------------------------------------------------

library(dplyr)
library(ggplot2)

# 1. Define all ventilation scenarios
ach_scenarios <- list(
  "All_0.5"      = c(0.5, 0.5, 0.5),
  "Seed6"        = c(6, 0.5, 0.5),
  "Other6"       = c(0.5, 6, 0.5),
  "Both_Other6"  = c(0.5, 6, 6),
  "Seed+Other6"  = c(6, 6, 0.5),
  "All_6"        = c(6, 6, 6)
)

# 2. Run simulations for each scenario; aggregate total outbreak size per run
set.seed(123)
all_results <- lapply(names(ach_scenarios), function(scenario) {
  ach <- ach_scenarios[[scenario]]
  sims <- run_multiple_patch_simulations(num_runs = 10000, ACH_vals = ach)
  sims %>%
    group_by(run_id) %>%
    summarize(
      scenario = scenario,
      total_cases = sum(final_size),
      .groups = "drop"
    )
})
final_df <- bind_rows(all_results)

# 3. Build scenario-level summary columns (mean ACH, was seed patch high ACH?)
scenario_design <- data.frame(
  scenario = names(ach_scenarios),
  mean_ACH = sapply(ach_scenarios, mean),
  seed_6   = sapply(ach_scenarios, function(x) as.integer(x[1] == 6)),
  other6   = sapply(ach_scenarios, function(x) sum(x[-1] == 6))
)

# 4. Merge scenario design info into results
final_df <- final_df %>% left_join(scenario_design, by = "scenario")

# ---- EXTINCTION ANALYSIS ----
# Set extinction threshold (change as needed: 1=stringent, 5/10/20=relaxed)
extinction_threshold <- 10   # Set to 5, 10, 20 to try alternatives

# 5. Summarize scenario-level stats, including extinction probability
summary_stats <- final_df %>%
  group_by(scenario, mean_ACH, seed_6) %>%
  summarize(
    mean_size = mean(total_cases),
    median_size = median(total_cases),
    ci_lower = quantile(total_cases, 0.025),
    ci_upper = quantile(total_cases, 0.975),
    extinction_rate = mean(total_cases <= extinction_threshold),
    .groups = "drop"
  )
print(summary_stats)

# ---- PLOTTING ----

# 6. Violin+boxplot of total outbreak size by scenario (color by seed_6)
ggplot(final_df, aes(x = scenario, y = total_cases, fill = factor(seed_6))) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("0" = "salmon", "1" = "skyblue"),
                    labels = c("Seed Patch <6", "Seed Patch 6"),
                    name = "Seed Patch\nVentilation") +
  labs(
    title = "Outbreak Sizes by ACH Scenario (Colored by Seed Patch Ventilation)",
    x = "Scenario",
    y = "Total Cases (All Patches)"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 7. Mean outbreak size vs mean ACH, colored by seed_6 (NO error bars)
ggplot(summary_stats, aes(x = mean_ACH, y = mean_size, color = factor(seed_6))) +
  geom_point(size = 3) +
  labs(
    title = "Mean Outbreak Size vs Mean ACH\n(Colored by Seed Patch Ventilation)",
    x = "Mean Air Changes per Hour (ACH)",
    y = "Mean Outbreak Size",
    color = "Seed Patch\nVentilation"
  ) +
  scale_color_manual(values = c("0" = "salmon", "1" = "skyblue"),
                     labels = c("Seed Patch <6", "Seed Patch 6")) +
  theme_minimal(base_size = 14)

# 8. Extinction probability by scenario (with your chosen threshold)
library(scales)
ggplot(summary_stats, aes(x = reorder(scenario, extinction_rate), y = extinction_rate, color = scenario)) +
  geom_point(size = 3) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = paste("Extinction Probability by Scenario\n(Threshold <=", extinction_threshold, "cases)"),
    x = "Scenario",
    y = "Extinction Rate"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# 9. Mean outbreak size by scenario (no error bars, by scenario label)
ggplot(summary_stats, aes(x = reorder(scenario, mean_size), y = mean_size, color = scenario)) +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    title = "Mean Outbreak Size by Scenario",
    x = "Scenario",
    y = "Mean Outbreak Size"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# -------------------------------------------------
# BLOCK 5: Regression Models & Diagnostics
# -------------------------------------------------
# Fit both negative binomial (for outbreak size) and logistic (for extinction) models
# using mean ACH and seed patch ventilation as predictors.
#
# Assumes 'final_df' from Block 4 (with scenario design columns merged).

library(MASS)    # for glm.nb
library(dplyr)
library(ggplot2)
library(scales)

# 1. (If not already present) Build scenario-level design info
scenario_design <- data.frame(
  scenario = names(ach_scenarios),
  mean_ACH = sapply(ach_scenarios, mean),
  seed_6   = sapply(ach_scenarios, function(x) as.integer(x[1] == 6)),
  other6   = sapply(ach_scenarios, function(x) sum(x[-1] == 6))
)

# 2. Merge predictors, calculate outcomes
model_df <- final_df %>%
  left_join(scenario_design, by = "scenario") %>%
  mutate(
    adj_size = pmax(0, total_cases - 1),           # Secondary cases (not <0)
    extinct = as.integer(total_cases <= 1)          # 1 = extinction, 0 = outbreak
  )

# 3. Negative binomial regression on outbreak size (secondary cases)
model_nb <- glm.nb(adj_size ~ mean_ACH.x + seed_6.x, data = model_df)
cat("\n--- Negative Binomial Regression ---\n")
print(summary(model_nb))

# 4. Logistic regression for extinction probability
logit_mod <- glm(extinct ~ mean_ACH.x + seed_6.x, data = model_df, family = binomial())
cat("\n--- Logistic Regression for Extinction ---\n")
print(summary(logit_mod))

# 5. Predict and plot secondary outbreak size (NegBin)
pred_data_nb <- expand.grid(
  mean_ACH.x = seq(0.5, 6, by = 0.1),
  seed_6.x   = c(0, 1)
)
pred_data_nb$predicted <- predict(model_nb, newdata = pred_data_nb, type = "response")

ggplot(pred_data_nb, aes(x = mean_ACH.x, y = predicted, color = factor(seed_6.x))) +
  geom_line(size = 1.2) +
  labs(
    title = "Predicted Secondary Outbreak Size (NegBin Model)",
    x = "Mean Air Changes per Hour (ACH)",
    y = "Predicted Secondary Cases",
    color = "Seed Patch High ACH"
  ) +
  scale_color_manual(values = c("red", "blue"), labels = c("No", "Yes")) +
  theme_minimal(base_size = 14)

# 6. Predict and plot extinction probability (logistic)
pred_data_logit <- expand.grid(
  mean_ACH.x = seq(0.5, 6, by = 0.1),
  seed_6.x   = c(0, 1)
)
pred_data_logit$prob_extinct <- predict(logit_mod, newdata = pred_data_logit, type = "response")

ggplot(pred_data_logit, aes(x = mean_ACH.x, y = prob_extinct, color = factor(seed_6.x))) +
  geom_line(size = 1.2) +
  labs(
    title = "Predicted Extinction Probability (Logistic Model)",
    x = "Mean Air Changes per Hour (ACH)",
    y = "Probability of Extinction (≤10 case)",
    color = "Seed Patch High ACH"
  ) +
  scale_color_manual(values = c("red", "blue"), labels = c("No", "Yes")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 14)

# 7. (Optional) Diagnostics for NegBin model
res_nb <- residuals(model_nb, type = "pearson")
fit_nb <- fitted(model_nb)
plot(fit_nb, res_nb, xlab = "Fitted", ylab = "Pearson Residuals",
     main = "NegBin Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)
qqnorm(res_nb)
qqline(res_nb, col = "red")
hist(res_nb, breaks = 30, col = "skyblue",
     main = "Histogram of Residuals", xlab = "Pearson Residuals")

# -------------------------------------------------
# BLOCK 6: Parameter Sweep on "Good" Ventilation ACH
# -------------------------------------------------

good_ACH_seq <- seq(2, 12, by = 1)   # Or use by = 0.5 for finer granularity

# Storage for results across sweeps
sweep_summary_list <- list()
sweep_model_list <- list()

for (good_ACH in good_ACH_seq) {
  # Redefine scenarios, replacing '6' with current 'good_ACH'
  ach_scenarios_sweep <- list(
    "All_Low"      = rep(0.5, 3),
    "SeedGood"     = c(good_ACH, 0.5, 0.5),
    "OtherGood"    = c(0.5, good_ACH, 0.5),
    "BothOtherGood"= c(0.5, good_ACH, good_ACH),
    "Seed+Other"   = c(good_ACH, good_ACH, 0.5),
    "All_Good"     = rep(good_ACH, 3)
  )
  
  # Run simulations as before
  set.seed(123)
  all_results <- lapply(names(ach_scenarios_sweep), function(scenario) {
    ach <- ach_scenarios_sweep[[scenario]]
    sims <- run_multiple_patch_simulations(num_runs = 10000, ACH_vals = ach)
    sims %>%
      group_by(run_id) %>%
      summarize(
        scenario = scenario,
        total_cases = sum(final_size),
        .groups = "drop"
      )
  })
  final_df <- bind_rows(all_results)
  
  # Add design columns
  scenario_design <- data.frame(
    scenario = names(ach_scenarios_sweep),
    mean_ACH = sapply(ach_scenarios_sweep, mean),
    seed_good = sapply(ach_scenarios_sweep, function(x) as.integer(x[1] == good_ACH)),
    other_good = sapply(ach_scenarios_sweep, function(x) sum(x[-1] == good_ACH))
  )
  
  final_df <- final_df %>% left_join(scenario_design, by = "scenario")
  
  # Summarize outcomes for each scenario at this ACH
  extinction_threshold <- 1
  summary_stats <- final_df %>%
    group_by(scenario, mean_ACH, seed_good) %>%
    summarize(
      mean_size = mean(total_cases),
      median_size = median(total_cases),
      extinction_rate = mean(total_cases <= extinction_threshold),
      .groups = "drop"
    )
  summary_stats$good_ACH <- good_ACH  # Track the current good ACH
  
  # Fit NB model (on all runs pooled across scenarios)
  model_df <- final_df %>%
    mutate(adj_size = pmax(0, total_cases - 1))
  model_nb <- MASS::glm.nb(adj_size ~ mean_ACH + seed_good, data = model_df)
  
  # Store
  sweep_summary_list[[as.character(good_ACH)]] <- summary_stats
  sweep_model_list[[as.character(good_ACH)]] <- coef(model_nb)
}

# Combine sweep summaries
sweep_summary <- dplyr::bind_rows(sweep_summary_list, .id = "good_ACH")
sweep_summary$good_ACH <- as.numeric(sweep_summary$good_ACH)

# Extract effect size for "seed patch = good ACH" from NB models
nb_effects <- data.frame(
  good_ACH = as.numeric(names(sweep_model_list)),
  seed_good_effect = sapply(sweep_model_list, function(coefs) coefs["seed_good"])
)

# --- Plot: NB effect of seed patch ventilation vs 'good' ACH threshold ---
library(ggplot2)
ggplot(nb_effects, aes(x = good_ACH, y = seed_good_effect)) +
  geom_line(size = 1.3, color = "blue") +
  geom_point(size = 2, color = "blue") +
  labs(
    title = "Effect of 'Good' Ventilation in Seed Patch vs. Definition of Good (ACH)",
    x = "Definition of 'Good' Ventilation (ACH)",
    y = "NB Regression Coefficient for Seed Patch = Good"
  ) +
  theme_minimal(base_size = 15)

# --- Plot: Mean outbreak size vs mean ACH, stratified by seed patch ---
ggplot(sweep_summary, aes(x = mean_ACH, y = mean_size, color = factor(seed_good))) +
  geom_point(size = 2.5) +
  geom_line(aes(group = interaction(seed_good, good_ACH)), alpha = 0.5) +
  facet_wrap(~ good_ACH, ncol = 2, scales = "free_x") +
  labs(
    title = "Mean Outbreak Size vs. Mean ACH (across Good ACH definitions)",
    x = "Mean ACH in Scenario",
    y = "Mean Outbreak Size",
    color = "Seed Patch Good?"
  ) +
  theme_minimal(base_size = 13)

# --- Plot: Extinction probability vs mean ACH, stratified by seed patch ---
ggplot(sweep_summary, aes(x = mean_ACH, y = extinction_rate, color = factor(seed_good))) +
  geom_point(size = 2.5) +
  geom_line(aes(group = interaction(seed_good, good_ACH)), alpha = 0.5) +
  facet_wrap(~ good_ACH, ncol = 2, scales = "free_x") +
  labs(
    title = "Extinction Rate vs. Mean ACH (across Good ACH definitions)",
    x = "Mean ACH in Scenario",
    y = "Extinction Probability",
    color = "Seed Patch Good?"
  ) +
  theme_minimal(base_size = 13)



############################################
# VAPOR: SINGLE-COMPARTMENT DETERMINISTIC MODEL
# Includes joint and de-overlapped transmission logic,
# route attribution, and R₀ sensitivity analysis
# Author: David Fisman et al.
# Updated: June 2025
############################################

# ---------------------------
# Core Parameters
# ---------------------------
N <- 100                   # Total population
q1 <- 0.005                # Close-contact transmission probability
q <- 50                    # Quanta generation rate (quanta/hour)
p <- 0.5                   # Pulmonary ventilation rate (m3/hour)
V <- 100                   # Room volume (m3)
ACH <- 2                   # Air changes per hour
dt <- 1                    # Duration of one generation (hour)
initial_infected <- 1      # Initial cases
max_time <- 50             # Simulation duration (generations)
f <- 0.1                   # Fraction of aerosolizing cases (Cb)

# ---------------------------
# Option 1: Combined Probability Transmission
# ---------------------------
run_combined_model <- function(max_time = 50) {
  S <- N - initial_infected
  C <- initial_infected
  incidence <- numeric(max_time + 1)
  incidence[1] <- C
  
  for (t in 1:max_time) {
    if (S == 0 || C == 0) break
    Cb <- max(1, round(f * C))
    Ca <- max(0, C - Cb)
    
    P_c <- 1 - (1 - q1)^(Ca + Cb)
    P_a <- 1 - exp(-Cb * q * p * dt / (V * ACH))
    P_total <- 1 - (1 - P_c) * (1 - P_a)
    
    new_C <- min(S, round(S * P_total))
    S <- S - new_C
    C <- new_C
    incidence[t + 1] <- C
  }
  return(incidence)
}

# ---------------------------
# Option 2: De-overlapped with Route Attribution
# ---------------------------
run_deoverlapped_with_attribution <- function(max_time = 50) {
  S <- N - initial_infected
  C <- initial_infected
  incidence <- numeric(max_time + 1)
  incidence_contact <- numeric(max_time + 1)
  incidence_airborne <- numeric(max_time + 1)
  
  incidence[1] <- C
  
  for (t in 1:max_time) {
    if (S == 0 || C == 0) break
    Cb <- max(1, round(f * C))
    Ca <- max(0, C - Cb)
    
    P_c <- 1 - (1 - q1)^(Ca + Cb)
    P_a <- 1 - exp(-Cb * q * p * dt / (V * ACH))
    
    E_c <- S * P_c
    E_a <- (S - E_c) * P_a
    
    new_Cc <- round(E_c)
    new_Ca <- round(E_a)
    new_C <- min(S, new_Cc + new_Ca)
    
    # Proportional rescaling to prevent overflow
    if ((new_Cc + new_Ca) > S) {
      scale <- S / (new_Cc + new_Ca)
      new_Cc <- round(new_Cc * scale)
      new_Ca <- round(new_Ca * scale)
      new_C <- new_Cc + new_Ca
    }
    
    S <- S - new_C
    C <- new_C
    incidence[t + 1] <- new_C
    incidence_contact[t + 1] <- new_Cc
    incidence_airborne[t + 1] <- new_Ca
  }
  
  return(list(
    total = incidence,
    contact = incidence_contact,
    airborne = incidence_airborne
  ))
}

# ---------------------------
# Plot Attribution Results
# ---------------------------
result <- run_deoverlapped_with_attribution()
barplot(rbind(result$contact, result$airborne),
        beside = FALSE,
        col = c("orange", "skyblue"),
        names.arg = 0:max_time,
        xlab = "Generation",
        ylab = "New Cases",
        main = "Route Attribution: Contact vs Airborne",
        border = NA,
        legend.text = c("Contact", "Airborne"),
        args.legend = list(x = "topright", bty = "n"))

# ---------------------------
# Analytical R₀ Calculation
# ---------------------------
compute_R0 <- function(q1, q, p, dt, V, ACH, f = 0.1, N = 100) {
  P_WR <- 1 - exp(-q * p * dt / (V * ACH))
  R0_Ca <- (1 - f) * N * q1
  R0_Cb <- f * N * (q1 + P_WR)
  return(R0_Ca + R0_Cb)
}

# ---------------------------
# Sensitivity Analyses: R₀ vs. ACH
# ---------------------------
ACH_vals <- seq(0.5, 14, by = 0.5)

# Sensitivity: Aerosolizing fraction f
f_vals <- c(0.05, 0.1, 0.2)
colors_f <- c("blue", "green", "purple")
plot(NULL, xlim = c(0.5, 14), ylim = c(0, 5),
     xlab = "ACH", ylab = expression(R[0]),
     main = "R₀ vs. ACH for Different Aerosol Fractions")
for (i in seq_along(f_vals)) {
  lines(ACH_vals,
        sapply(ACH_vals, function(x) compute_R0(q1, q, p, dt, V, x, f = f_vals[i])),
        col = colors_f[i], lwd = 2)
}
legend("topright", legend = paste("f =", f_vals), col = colors_f, lty = 1, lwd = 2)
abline(h = 1, col = "red", lty = 2)

# Sensitivity: Quanta generation rate q
q_vals <- c(25, 50, 100)
colors_q <- c("darkorange", "darkgreen", "darkred")
plot(NULL, xlim = c(0.5, 14), ylim = c(0, 8),
     xlab = "ACH", ylab = expression(R[0]),
     main = "R₀ vs. ACH for Different Room Volumes")
for (i in seq_along(q_vals)) {
  lines(ACH_vals,
        sapply(ACH_vals, function(x) compute_R0(q1, q_vals[i], p, dt, V, x, f = 0.1)),
        col = colors_q[i], lwd = 2)
}
legend("topright", legend = paste("q =", q_vals), col = colors_q, lty = 1, lwd = 2)
abline(h = 1, col = "red", lty = 2)

# Sensitivity: Room volume V
V_vals <- c(100, 250, 500)
colors_V <- c("brown", "navy", "darkcyan")
plot(NULL, xlim = c(0.5, 14), ylim = c(0, 6),
     xlab = "ACH", ylab = expression(R[0]),
     main = "R0 vs. ACH for Different Room Volumes")
for (i in seq_along(V_vals)) {
  lines(ACH_vals,
        sapply(ACH_vals, function(x) compute_R0(q1, q, p, dt, V_vals[i], x, f = 0.1)),
        col = colors_V[i], lwd = 2)
}
legend("topright", legend = paste("V =", V_vals), col = colors_V, lty = 1, lwd = 2)
abline(h = 1, col = "red", lty = 2)

#----------------------------------------------------------------------------------
#----------------------------MAKE ABOVE FIGURES WITH GGPLOT------------------------
#----------------------------------------------------------------------------------
library(ggplot2)
# Create a tidy data frame for ggplot
df_f <- do.call(rbind, lapply(f_vals, function(fv) {
  data.frame(
    ACH = ACH_vals,
    R0 = sapply(ACH_vals, function(x) compute_R0(q1, q, p, dt, V, x, f = fv)),
    f = factor(fv)
  )
}))

ggplot(df_f, aes(x = ACH, y = R0, color = f)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = expression(R[0] ~ "vs. ACH for Different Aerosol Fractions"),
    x = "Air Changes per Hour (ACH)",
    y = expression(R[0]),
    color = "Aerosol fraction (f)"
  ) +
  theme_minimal()

df_q <- do.call(rbind, lapply(q_vals, function(qv) {
  data.frame(
    ACH = ACH_vals,
    R0 = sapply(ACH_vals, function(x) compute_R0(q1, qv, p, dt, V, x, f = 0.1)),
    q = factor(qv)
  )
}))

# ---------------------------
# Sensitivity Plot: R₀ vs. ACH for Different Quanta Generation Rates
# ---------------------------

library(ggplot2)

# Prepare data
df_q <- do.call(rbind, lapply(q_vals, function(qv) {
  data.frame(
    ACH = ACH_vals,
    R0 = sapply(ACH_vals, function(x) compute_R0(q1, qv, p, dt, V, x, f = 0.1)),
    q = factor(qv)
  )
}))

# Create plot
p_q <- ggplot(df_q, aes(x = ACH, y = R0, color = q)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    #title = expression(R[0] ~ "vs. ACH for Different Quanta Generation Rates"),
    x = "Air Changes per Hour (ACH)",
    y = expression(R[0]),
    color = "Quanta Rate (q)"
  ) +
  theme_minimal(base_size = 14)

# Display plot
print(p_q)

# ---------------------------
# Save High-Quality Figures
# ---------------------------

# File names
pdf_file <- "Figure1_R0_vs_ACH_q_sensitivity.pdf"
png_file <- "Figure1_R0_vs_ACH_q_sensitivity.png"

# Save as high-resolution PDF
ggsave(filename = pdf_file, plot = p_q, width = 7, height = 5, units = "in", dpi = 300)

# Save as high-resolution PNG
ggsave(filename = png_file, plot = p_q, width = 7, height = 5, units = "in", dpi = 300)


df_V <- do.call(rbind, lapply(V_vals, function(Vv) {
  data.frame(
    ACH = ACH_vals,
    R0 = sapply(ACH_vals, function(x) compute_R0(q1, q, p, dt, Vv, x, f = 0.1)),
    V = factor(Vv)
  )
}))

ggplot(df_V, aes(x = ACH, y = R0, color = V)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = expression(R[0] ~ "vs. ACH for Different Room Volumes"),
    x = "Air Changes per Hour (ACH)",
    y = expression(R[0]),
    color = "Room Volume (m³)"
  ) +
  theme_minimal()

result <- run_deoverlapped_with_attribution()

df_attr <- data.frame(
  Generation = rep(0:max_time, 2),
  Cases = c(result$contact, result$airborne),
  Route = rep(c("Contact", "Airborne"), each = max_time + 1)
)

library(ggplot2)

# Build the stacked bar plot
p_route <- ggplot(df_attr, aes(x = Generation, y = Cases, fill = Route)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Contact" = "orange", "Airborne" = "skyblue")) +
  scale_x_continuous(limits = c(0, 15), expand = expansion(mult = c(0, 0.02))) +
  labs(
    #title = "Figure 1b: Route Attribution – Contact vs. Airborne Transmission",
    x = "Generation",
    y = "New Cases",
    fill = "Transmission Route"
  ) +
  theme_minimal(base_size = 16) +  # this sets base font size
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 15)
  )

# Display it
print(p_route)

# Export as high-resolution PDF and PNG
ggsave("Figure1b_RouteAttribution_Barplot.pdf", plot = p_route, width = 7, height = 5, units = "in", dpi = 300)
ggsave("Figure1b_RouteAttribution_Barplot.png", plot = p_route, width = 7, height = 5, units = "in", dpi = 300)

#-------------------------
#Create a composite figure
#-------------------------

# Rebuild df_q correctly
df_q <- do.call(rbind, lapply(q_vals, function(qv) {
  data.frame(
    ACH = ACH_vals,
    R0 = sapply(ACH_vals, function(x) compute_R0(q1, qv, p, dt, V, x, f = 0.1)),
    q = factor(qv)
  )
}))

# Redefine p_q without referencing ParamValue
p_q <- ggplot(df_q, aes(x = ACH, y = R0, color = q)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "Air Changes per Hour (ACH)",
    y = expression(R[0]),
    color = "Quanta Rate (q)"
  ) +
  theme_minimal(base_size = 16)


library(ggplot2)
library(patchwork)

# Panel a (ACH vs R₀, varying quanta)
p1a <- p_q +
  labs(title = NULL) +  # remove internal title
  theme(
    plot.margin = margin(10, 10, 10, 10),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

# Panel b (route attribution)
p1b <- p_route +
  labs(title = NULL) +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  )

# Compose vertically with labels
fig1 <- (p1a / p1b) +
  plot_annotation(
    tag_levels = 'a',
    #title = "Figure 1 | Impact of ventilation and route attribution on outbreak dynamics",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0),
      plot.tag = element_text(size = 16, face = "bold")
    )
  )

print(fig1)

# Save high-quality composite figure
ggsave("Figure1_Composite_VAPOR.pdf", plot = fig1, width = 7.5, height = 10, units = "in", dpi = 300)
ggsave("Figure1_Composite_VAPOR.png", plot = fig1, width = 7.5, height = 10, units = "in", dpi = 300)


#------------------------------------------------------------------------------------------
#Other sensitivity analyses for supplementary appendix-------------------------------------
#------------------------------------------------------------------------------------------

library(patchwork)

# p_f (ACH vs. R0 by aerosol fraction f)
p_f <- ggplot(df_f, aes(x = ACH, y = R0, color = f)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "ACH", y = expression(R[0]), color = "Aerosol fraction (f)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(10, 10, 10, 10))

# p_V (ACH vs. R0 by room volume V)
p_V <- ggplot(df_V, aes(x = ACH, y = R0, color = V)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "ACH", y = NULL, color = "Room volume (m³)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(10, 10, 10, 10))

# Combine side-by-side
fig_S1 <- (p_f | p_V) +
  plot_annotation(
    title = "Supplementary Figure S1 | Sensitivity of R₀ to airborne transmission parameters",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      plot.tag = element_text(size = 14, face = "bold")
    )
  )

print(fig_S1)
# Save
ggsave("FigureS1_Sensitivity_f_V.pdf", plot = fig_S1, width = 10, height = 4.5, units = "in", dpi = 300)
ggsave("FigureS1_Sensitivity_f_V.png", plot = fig_S1, width = 10, height = 4.5, units = "in", dpi = 300)

#library(patchwork)

# p_f (ACH vs. R0 by aerosol fraction f)
p_f <- ggplot(df_f, aes(x = ACH, y = R0, color = f)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "ACH", y = expression(R[0]), color = "Aerosol fraction (f)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(10, 10, 10, 10))

# p_V (ACH vs. R0 by room volume V)
p_V <- ggplot(df_V, aes(x = ACH, y = R0, color = V)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "ACH", y = expression(R[0]), color = "Room volume (m³)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(10, 10, 10, 10))

# Combine side-by-side
fig_S1 <- (p_f / p_V) +
  plot_annotation(
    #title = expression("Supplementary Figure S1 | Sensitivity of " * R[0] * " to airborne transmission parameters"),
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      plot.tag = element_text(size = 14, face = "bold")
    )
  )

print(fig_S1)
# Save
ggsave("FigureS1_Sensitivity_f_V.pdf", plot = fig_S1, width = 10, height = 4.5, units = "in", dpi = 300)
ggsave("FigureS1_Sensitivity_f_V.png", plot = fig_S1, width = 10, height = 4.5, units = "in", dpi = 300)


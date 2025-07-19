library(ggplot2)
library(reshape2)

# 1. Efron's Biased Coin Design (fixed p)
efron_bcd <- function(Dn, p = 2/3) {
  if (Dn < 0) return(p)
  if (Dn == 0) return(0.5)
  return(1 - p)
}

# 2. GBCD (Generalized Biased Coin Design)
gbcd_probability <- function(Dn, n, zeta = 1) {
  x <- Dn / n
  x <- max(min(x, 1), -1)  # Clip x to [-1,1]
  numerator <- (1 - x)^zeta
  denominator <- numerator + (1 + x)^zeta
  prob <- numerator / denominator
  return(prob)
}

# 3. ABCD (Adjustable Biased Coin Design)
abcd_f <- function(Dn, a = 1) {
  if (abs(Dn) <= 1) {
    return(0.5)
  } else if (Dn > 1) {
    return(1 / (1 + (Dn^a)))
  } else { 
    return((abs(Dn)^a) / (1 + (abs(Dn)^a)))
  }
}

# 4. Bayesian BCD with clipping and epsilon to avoid division by zero
bayesian_bcd <- function(Dn, n, gamma) {
  epsilon <- 1e-8
  x <- Dn / n
  x <- max(min(x, 1 - epsilon), -1 + epsilon)
  
  base_term <- 1 + (1 - x) / (n * (1 + x + epsilon))
  alt_term <- 1 + (1 + x) / (n * (1 - x + epsilon))
  
  numerator <- base_term^(1 / gamma)
  denominator <- numerator + alt_term^(1 / gamma)
  prob <- numerator / denominator
  return(prob)
}

# 5. Dominant BCD
dbcd_probability <- function(Dn, n) {
  a_n <- 1 / (floor(n / 20) + 1)
  if (abs(Dn) <= 1) {
    prob <- 0.5
  } else if (Dn > 1) {
    prob <- 1 / (1 + (Dn^a_n))
  } else { 
    prob <- (abs(Dn)^a_n) / (1 + (abs(Dn)^a_n))
  }
  return(prob)
}

# Prepare sequences
Dn_values <- seq(-10, 10, length.out = 300)
Dn_values_norm <- seq(-1, 1, length.out = 300)
n_fixed <- 50

# --- 1. GBCD plot ---
df_gbcd <- data.frame(
  Dn_norm = Dn_values_norm,
  GBCD_zeta1 = sapply(Dn_values_norm * n_fixed, function(d) gbcd_probability(d, n_fixed, zeta = 1)),
  GBCD_zeta2 = sapply(Dn_values_norm * n_fixed, function(d) gbcd_probability(d, n_fixed, zeta = 2)),
  GBCD_zeta4 = sapply(Dn_values_norm * n_fixed, function(d) gbcd_probability(d, n_fixed, zeta = 4))
)

df_gbcd_melt <- melt(df_gbcd, id.vars = "Dn_norm")

plot_gbcd <- ggplot(df_gbcd_melt, aes(x = Dn_norm, y = value, color = variable)) +
  geom_line(size = 1.1) +
  labs(title = "Generalized Biased Coin Design (GBCD)",
       x = "Normalized Imbalance (Dn/n)",
       y = "Allocation Probability",
       color = "Zeta") +
  theme_minimal() +
  theme(legend.position = "bottom")

# --- 2. Bayesian BCD plot ---
df_bayesian <- data.frame(
  Dn_norm = Dn_values_norm,
  Bayesian_gamma1 = sapply(Dn_values_norm * n_fixed, function(d) bayesian_bcd(d, n_fixed, gamma = 1)),
  Bayesian_gamma0.1 = sapply(Dn_values_norm * n_fixed, function(d) bayesian_bcd(d, n_fixed, gamma = 0.1)),
  Bayesian_gamma0.01 = sapply(Dn_values_norm * n_fixed, function(d) bayesian_bcd(d, n_fixed, gamma = 0.01))
)

df_bayesian_melt <- melt(df_bayesian, id.vars = "Dn_norm")

plot_bayesian <- ggplot(df_bayesian_melt, aes(x = Dn_norm, y = value, color = variable)) +
  geom_line(size = 1.1) +
  labs(title = "Bayesian Biased Coin Design (BBCD)",
       x = "Normalized Imbalance (Dn/n)",
       y = "Allocation Probability",
       color = "Gamma") +
  theme_minimal() +
  theme(legend.position = "bottom")

# --- 3. Efron vs ABCD plot ---
abcd_as <- c(0.5, 1, 2, 5)
df_abcd <- data.frame(Dn = Dn_values,
                      Efron_BCD = sapply(Dn_values, efron_bcd))

for (a_val in abcd_as) {
  col_name <- paste0("ABCD_a", a_val)
  df_abcd[[col_name]] <- sapply(Dn_values, function(d) abcd_f(d, a = a_val))
}

df_abcd_melt <- melt(df_abcd, id.vars = "Dn")

plot_abcd <- ggplot(df_abcd_melt, aes(x = Dn, y = value, color = variable)) +
  geom_line(size = 1.1) +
  labs(title = "Efron's BCD vs Adjustable BCD (ABCD)",
       x = "Imbalance (Dn)",
       y = "Allocation Probability",
       color = "Design") +
  theme_minimal() +
  theme(legend.position = "bottom")

# --- 4. Dominant BCD for different n ---
dominant_ns <- c(20, 40, 60, 80, 100)
df_dbcd <- data.frame(Dn = Dn_values)

for (n_val in dominant_ns) {
  col_name <- paste0("DBCD_n", n_val)
  df_dbcd[[col_name]] <- sapply(Dn_values, function(d) dbcd_probability(d, n_val))
}

df_dbcd_melt <- melt(df_dbcd, id.vars = "Dn")

plot_dbcd <- ggplot(df_dbcd_melt, aes(x = Dn, y = value, color = variable)) +
  geom_line(size = 1.1) +
  labs(title = "Dominant BCD Allocation Probabilities\nfor Different Sample Sizes (n)",
       x = "Imbalance (Dn)",
       y = "Allocation Probability",
       color = "Sample Size (n)") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Print all plots
print(plot_gbcd)
print(plot_bayesian)
print(plot_abcd)
print(plot_dbcd)

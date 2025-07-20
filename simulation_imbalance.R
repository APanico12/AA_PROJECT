library(ggplot2)
library(gridExtra)
library(grid)

# Allocation functions
efron_bcd <- function(Dn, p = 2/3) {
  if (Dn < 0) return(p)
  if (Dn == 0) return(0.5)
  return(1 - p)
}

gbcd_probability <- function(Dn, n, zeta = 1) {
  x <- Dn / n
  x <- max(min(x, 1), -1)
  numerator <- (1 - x)^zeta
  denominator <- numerator + (1 + x)^zeta
  numerator / denominator
}

abcd_f <- function(Dn, a = 1) {
  if (abs(Dn) <= 1) {
    return(0.5)
  } else if (Dn > 1) {
    return(1 / (1 + (Dn^a)))
  } else { 
    return((abs(Dn)^a) / (1 + (abs(Dn)^a)))
  }
}

bayesian_bcd <- function(Dn, n, gamma) {
  epsilon <- 1e-8
  x <- Dn / n
  x <- max(min(x, 1 - epsilon), -1 + epsilon)
  base_term <- 1 + (1 - x) / (n * (1 + x + epsilon))
  alt_term <- 1 + (1 + x) / (n * (1 - x + epsilon))
  numerator <- base_term^(1 / gamma)
  denominator <- numerator + alt_term^(1 / gamma)
  numerator / denominator
}

dbcd_probability <- function(Dn, n) {
  a_n <- 1 / (floor(n / 20) + 1)
  if (abs(Dn) <= 1) {
    0.5
  } else if (Dn > 1) {
    1 / (1 + (Dn^a_n))
  } else {
    (abs(Dn)^a_n) / (1 + (abs(Dn)^a_n))
  }
}

# Designs list with parameters
designs <- list(
  list(name = "BCD (2/3)", fun = efron_bcd, params = list(p = 2/3)),
  list(name = "ABCD[F(1)]", fun = abcd_f, params = list(a = 1)),
  list(name = "GBCD (ζ=1)", fun = gbcd_probability, params = list(zeta = 1)),
  list(name = "GBCD (ζ=2)", fun = gbcd_probability, params = list(zeta = 2)),
  list(name = "BBCD (γ=1)", fun = bayesian_bcd, params = list(gamma = 1)),
  list(name = "BBCD (γ=0.1)", fun = bayesian_bcd, params = list(gamma = 0.1)),
  list(name = "DomBCD", fun = dbcd_probability, params = list())
)

sample_sizes <- c(20, 40, 60, 80)

# Run simulation
results <- do.call(rbind, lapply(sample_sizes, function(n) {
  Dn_values <- seq(-n, n, by = 1)
  do.call(rbind, lapply(designs, function(design) {
    prob <- sapply(Dn_values, function(Dn) {
      if (design$name == "BCD (2/3)") {
        design$fun(Dn, p = design$params$p)
      } else if (design$name %in% c("GBCD (ζ=1)", "GBCD (ζ=2)")) {
        design$fun(Dn, n, zeta = design$params$zeta)
      } else if (design$name %in% c("BBCD (γ=1)", "BBCD (γ=0.1)")) {
        design$fun(Dn, n, gamma = design$params$gamma)
      } else if (design$name == "ABCD[F(1)]") {
        design$fun(Dn, a = design$params$a)
      } else if (design$name == "DomBCD") {
        design$fun(Dn, n)
      } else {
        NA
      }
    })
    data.frame(
      SampleSize = n,
      Design = design$name,
      Dn = Dn_values,
      Probability = prob
    )
  }))
}))

# Colors and linetypes for designs
colors <- c(
  "BCD (2/3)" = "red",
  "ABCD[F(1)]" = "blue",
  "GBCD (ζ=1)" = "green",
  "GBCD (ζ=2)" = "purple",
  "BBCD (γ=1)" = "orange",
  "BBCD (γ=0.1)" = "brown",
  "DomBCD" = "black"
)

linetypes <- c(
  "BCD (2/3)" = "solid",
  "ABCD[F(1)]" = "dashed",
  "GBCD (ζ=1)" = "dotted",
  "GBCD (ζ=2)" = "dotdash",
  "BBCD (γ=1)" = "twodash",
  "BBCD (γ=0.1)" = "longdash",
  "DomBCD" = "solid"
)

# Create plots without legends
plots_no_legend <- lapply(sample_sizes, function(n) {
  df <- subset(results, SampleSize == n)
  ggplot(df, aes(x = Dn, y = Probability, color = Design, linetype = Design)) +
    geom_line(size = 1) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    labs(title = paste("Sample Size n =", n),
         x = expression(D[n]),
         y = "Allocation Probability to Treatment A") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    ) +
    ylim(0, 1)
})

# Create a plot to extract legend
legend_plot <- ggplot(results, aes(x = Dn, y = Probability, color = Design, linetype = Design)) +
  geom_line() +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1.5, "lines"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))

# Function to extract legend grob
get_legend <- function(myplot) {
  tmp <- ggplot_gtable(ggplot_build(myplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

legend <- get_legend(legend_plot)

# Arrange plots in 2x2 grid
plots_grid <- arrangeGrob(grobs = plots_no_legend, nrow = 2, ncol = 2)

# Draw all together with legend on the right
grid.newpage()
grid.arrange(plots_grid, legend, ncol = 2, widths = c(4, 1))

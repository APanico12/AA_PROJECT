
# Adaptive Design Comparison – Full Simulation
rm(list = ls())

library(ggplot2)
library(ggpubr)
library(progress)

# Settings
n <- 200
M <- 5000
p <- 2/3
a <- c(3, 1, 1/3)
z.gbcd <- c(1, 2, 4)
gamma.bbayes <- c(1, 0.1)  

set.seed(1)

# Data containers
delta.Efron <- G.Efron <- D.Efron <- matrix(NA, ncol = n, nrow = length(p))
delta.ABCD <- G.ABCD <- D.ABCD <- matrix(NA, ncol = n, nrow = length(a))
delta.Dom <- G.Dom <- D.Dom <- matrix(NA, ncol = n, nrow = 1)
delta.GBCD <- G.GBCD <- D.GBCD <- matrix(NA, ncol = n, nrow = length(z.gbcd))
delta.BBayes2 <- G.BBayes2 <- D.BBayes2 <- matrix(NA, ncol = n, nrow = length(gamma.bbayes))


SB.Efron <- L.Efron <- array(NA, dim = c(M, n, length(p)))
SB.ABCD <- L.ABCD <- array(NA, dim = c(M, n, length(a)))
SB.Dom <- L.Dom <- array(NA, dim = c(M, n, 1))
SB.GBCD <- L.GBCD <- array(NA, dim = c(M, n, length(z.gbcd)))
SB.BBayes2 <- L.BBayes2 <- array(NA, dim = c(M, n, length(gamma.bbayes)))


pb <- progress_bar$new(format = " [:bar] :percent eta: :eta", total = M, width = 60)

for (h in 1:M) {
  # Initial patient
  delta.Efron[, 1] <- rbinom(length(p), 1, .5)
  delta.ABCD[, 1] <- rbinom(length(a), 1, .5)
  delta.Atkin[, 1] <- rbinom(length(t.atkin), 1, .5)
  delta.Dom[, 1] <- rbinom(1, 1, .5)
  
  G.Efron[, 1] <- G.ABCD[, 1] <- rbinom(length(p), 1, .5)
  G.Dom[1, 1] <- rbinom(1, 1, .5)
  
  D.Efron[, 1] <- 2 * delta.Efron[, 1] - 1
  D.ABCD[, 1] <- 2 * delta.ABCD[, 1] - 1
  D.Atkin[, 1] <- 2 * delta.Atkin[, 1] - 1
  D.Dom[1, 1] <- 2 * delta.Dom[1, 1] - 1
  
  for (i in 2:n) {
    for (e in 1:length(p)) {
      if (D.Efron[e, i - 1] == 0) {
        delta.Efron[e, i] <- rbinom(1, 1, .5)
        G.Efron[e, i] <- rbinom(1, 1, .5)
      } else if (D.Efron[e, i - 1] >= 1) {
        delta.Efron[e, i] <- rbinom(1, 1, 1 - p[e])
        G.Efron[e, i] <- as.numeric(delta.Efron[e, i] == 0)
      } else {
        delta.Efron[e, i] <- rbinom(1, 1, p[e])
        G.Efron[e, i] <- as.numeric(delta.Efron[e, i] == 1)
      }
      D.Efron[e, i] <- 2 * sum(delta.Efron[e, 1:i]) - i
    }
    
    for (j in 1:length(a)) {
      if (D.ABCD[j, i - 1] == 0) {
        delta.ABCD[j, i] <- rbinom(1, 1, .5)
        G.ABCD[j, i] <- rbinom(1, 1, .5)
      } else if (D.ABCD[j, i - 1] >= 1) {
        delta.ABCD[j, i] <- rbinom(1, 1, 1 / (D.ABCD[j, i - 1]^a[j] + 1))
        G.ABCD[j, i] <- as.numeric(delta.ABCD[j, i] == 0)
      } else {
        delta.ABCD[j, i] <- rbinom(1, 1, abs(D.ABCD[j, i - 1])^a[j] / (abs(D.ABCD[j, i - 1])^a[j] + 1))
        G.ABCD[j, i] <- as.numeric(delta.ABCD[j, i] == 1)
      }
      D.ABCD[j, i] <- 2 * sum(delta.ABCD[j, 1:i]) - i
    }
    
    for (t in 1:length(t.atkin)) {
      if (D.Atkin[t, i - 1] == 0) {
        delta.Atkin[t, i] <- rbinom(1, 1, .5)
        G.Atkin[t, i] <- rbinom(1, 1, .5)
      } else {
        imbalance <- sum(1 - delta.Atkin[t, 1:(i - 1)])^t.atkin[t]
        balance <- sum(delta.Atkin[t, 1:(i - 1)])^t.atkin[t]
        p.alloc <- imbalance / (balance + imbalance)
        delta.Atkin[t, i] <- rbinom(1, 1, p.alloc)
        G.Atkin[t, i] <- as.numeric((D.Atkin[t, i - 1] >= 1 && delta.Atkin[t, i] == 0) ||
                                      (D.Atkin[t, i - 1] <= -1 && delta.Atkin[t, i] == 1))
      }
      D.Atkin[t, i] <- 2 * sum(delta.Atkin[t, 1:i]) - i
    }
    
    d.prev <- D.Dom[1, i - 1]
    if (abs(d.prev) <= 1) {
      delta.Dom[1, i] <- rbinom(1, 1, 0.5)
      G.Dom[1, i] <- rbinom(1, 1, 0.5)
    } else if (d.prev > 1) {
      delta.Dom[1, i] <- 0
      G.Dom[1, i] <- 1
    } else {
      delta.Dom[1, i] <- 1
      G.Dom[1, i] <- 1
    }
    D.Dom[1, i] <- 2 * sum(delta.Dom[1, 1:i]) - i
  }
  
  for (e in 1:length(p)) {
    SB.Efron[h, , e] <- 2 * cumsum(G.Efron[e, ]) / 1:n - 1
    L.Efron[h, , e] <- D.Efron[e, ]^2 / 1:n
  }
  for (j in 1:length(a)) {
    SB.ABCD[h, , j] <- 2 * cumsum(G.ABCD[j, ]) / 1:n - 1
    L.ABCD[h, , j] <- D.ABCD[j, ]^2 / 1:n
  }
  for (g in 1:length(z.gbcd)) {
    if (i == 2) {
      delta.GBCD[g, 1] <- rbinom(1, 1, 0.5)
      G.GBCD[g, 1] <- rbinom(1, 1, 0.5)
      D.GBCD[g, 1] <- 2 * delta.GBCD[g, 1] - 1
    }
    
    d.prev <- D.GBCD[g, i - 1]
    x <- d.prev / (i - 1)
    fz <- ((1 - x)^z.gbcd[g]) / ((1 - x)^z.gbcd[g] + (1 + x)^z.gbcd[g])
    
    delta.GBCD[g, i] <- rbinom(1, 1, fz)
    G.GBCD[g, i] <- as.numeric((d.prev >= 1 && delta.GBCD[g, i] == 0) ||
                                 (d.prev <= -1 && delta.GBCD[g, i] == 1))
    D.GBCD[g, i] <- 2 * sum(delta.GBCD[g, 1:i]) - i
  }
    for (g in 1:length(z.gbcd)) {
    SB.GBCD[h, , g] <- 2 * cumsum(G.GBCD[g, ]) / 1:n - 1
    L.GBCD[h, , g] <- D.GBCD[g, ]^2 / 1:n
    }
  for (g in 1:length(gamma.bbayes)) {
    if (i == 2) {
      delta.BBayes2[g, 1] <- rbinom(1, 1, 0.5)
      G.BBayes2[g, 1] <- rbinom(1, 1, 0.5)
      D.BBayes2[g, 1] <- 2 * delta.BBayes2[g, 1] - 1
    }
    
    d.prev <- D.BBayes2[g, i - 1]
    x <- d.prev / (i - 1)
    gamma <- gamma.bbayes[g]
    
    if (!is.na(x) && is.finite(x) && abs(x) < 1 && x != 0) {
      denom1 <- (i - 1) * (1 + x)
      denom2 <- (i - 1) * (1 - x)
      
      if (denom1 != 0 && denom2 != 0) {
        num <- (1 + (1 - x) / denom1)^(1 / gamma)
        den <- num + (1 + (1 + x) / denom2)^(1 / gamma)
        f.gamma <- num / den
      } else {
        f.gamma <- 0.5
      }
    } else {
      f.gamma <- 0.5
    }
    
    delta.BBayes2[g, i] <- rbinom(1, 1, f.gamma)
    G.BBayes2[g, i] <- as.numeric((d.prev >= 1 && delta.BBayes2[g, i] == 0) ||
                                    (d.prev <= -1 && delta.BBayes2[g, i] == 1))
    D.BBayes2[g, i] <- 2 * sum(delta.BBayes2[g, 1:i]) - i
  }
  for (g in 1:length(gamma.bbayes)) {
    SB.BBayes2[h, , g] <- 2 * cumsum(G.BBayes2[g, ]) / 1:n - 1
    L.BBayes2[h, , g] <- D.BBayes2[g, ]^2 / 1:n
  }
  
  
  SB.Dom[h, , 1] <- 2 * cumsum(G.Dom[1, ]) / 1:n - 1
  L.Dom[h, , 1] <- D.Dom[1, ]^2 / 1:n
  
  pb$tick()
}

# Averages
Loss.Efron <- apply(L.Efron, c(2, 3), mean)
Sel.Bias.Efron <- apply(SB.Efron, c(2, 3), mean)
Loss.ABCD <- apply(L.ABCD, c(2, 3), mean)
Sel.Bias.ABCD <- apply(SB.ABCD, c(2, 3), mean)
Loss.Dom <- colMeans(L.Dom[, , 1])
Sel.Bias.Dom <- colMeans(SB.Dom[, , 1])
Loss.GBCD <- apply(L.GBCD, c(2,3), mean)
Sel.Bias.GBCD <- apply(SB.GBCD, c(2,3), mean)
Loss.BBayes<- apply(L.BBayes2, c(2,3), mean)
Sel.Bias.BBayes <- apply(SB.BBayes2, c(2,3), mean)

# Labels and plot
lab.s <- c(
  paste("Efron(", round(p, 2), ")", sep = ""),
  paste("ABCD(", round(a, 2), ")", sep = ""),
  paste("GBCD(", z.gbcd, ")", sep = ""),
  paste("Bayesian(", gamma.bbayes, ")", sep = ""),
  "Dominant"
)

dfL <- data.frame(
  x = rep(1:n, length(lab.s)),
  y = c(as.vector(Loss.Efron), as.vector(Loss.ABCD), as.vector(Loss.GBCD),
        as.vector(Loss.BBayes), Loss.Dom),
  Design = factor(rep(lab.s, each = n))
)

dfB <- data.frame(
  x = rep(1:n, length(lab.s)),
  y = c(as.vector(Sel.Bias.Efron), as.vector(Sel.Bias.ABCD), as.vector(Sel.Bias.GBCD),
        as.vector(Sel.Bias.BBayes), Sel.Bias.Dom),
  Design = factor(rep(lab.s, each = n))
)

plotL <- ggplot(dfL, aes(x = x, y = y)) +
  geom_line(aes(linetype = Design, color = Design), size = 1) +
  labs(x = "n", y = "Loss") +
  theme_bw() + 
  theme(legend.position = 'bottom') +
  coord_cartesian(xlim = c(15, n))

plotB <- ggplot(dfB, aes(x = x, y = y)) +
  geom_line(aes(linetype = Design, color = Design), size = 1) +
  labs(x = "n", y = "Selection Bias") +
  theme_bw() + 
  theme(legend.position = 'bottom') +
  coord_cartesian(xlim = c(15, n))

# Show plots
windows()
#ggarrange(plotL, plotB, nrow = 2, common.legend = TRUE, legend = "bottom")
print(plotL)
print(plotB)

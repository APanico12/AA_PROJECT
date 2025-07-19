
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
t.atkin <- 1:2
wei.t <- 1:2

set.seed(1)

# Data containers
delta.Efron <- G.Efron <- D.Efron <- matrix(NA, ncol = n, nrow = length(p))
delta.ABCD <- G.ABCD <- D.ABCD <- matrix(NA, ncol = n, nrow = length(a))
delta.Atkin <- G.Atkin <- D.Atkin <- matrix(NA, ncol = n, nrow = length(t.atkin))
delta.Wei <- G.Wei <- D.Wei <- matrix(NA, ncol = n, nrow = length(wei.t))
delta.Dom <- G.Dom <- D.Dom <- matrix(NA, ncol = n, nrow = 1)

SB.Efron <- L.Efron <- array(NA, dim = c(M, n, length(p)))
SB.ABCD <- L.ABCD <- array(NA, dim = c(M, n, length(a)))
SB.Atkin <- L.Atkin <- array(NA, dim = c(M, n, length(t.atkin)))
SB.Wei <- L.Wei <- array(NA, dim = c(M, n, length(wei.t)))
SB.Dom <- L.Dom <- array(NA, dim = c(M, n, 1))

pb <- progress_bar$new(format = " [:bar] :percent eta: :eta", total = M, width = 60)

for (h in 1:M) {
  # Initial patient
  delta.Efron[, 1] <- rbinom(length(p), 1, .5)
  delta.ABCD[, 1] <- rbinom(length(a), 1, .5)
  delta.Atkin[, 1] <- rbinom(length(t.atkin), 1, .5)
  delta.Wei[, 1] <- rbinom(length(wei.t), 1, .5)
  delta.Dom[, 1] <- rbinom(1, 1, .5)
  
  G.Efron[, 1] <- G.ABCD[, 1] <- G.Atkin[, 1] <- G.Wei[, 1] <- rbinom(length(p), 1, .5)
  G.Dom[1, 1] <- rbinom(1, 1, .5)
  
  D.Efron[, 1] <- 2 * delta.Efron[, 1] - 1
  D.ABCD[, 1] <- 2 * delta.ABCD[, 1] - 1
  D.Atkin[, 1] <- 2 * delta.Atkin[, 1] - 1
  D.Wei[, 1] <- 2 * delta.Wei[, 1] - 1
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
    
    for (w in 1:length(wei.t)) {
      x <- D.Wei[w, i - 1] / (i - 1)
      f.x <- (1 - x)^wei.t[w] / ((1 - x)^wei.t[w] + (1 + x)^wei.t[w])
      delta.Wei[w, i] <- rbinom(1, 1, f.x)
      G.Wei[w, i] <- as.numeric((D.Wei[w, i - 1] >= 1 && delta.Wei[w, i] == 0) ||
                                  (D.Wei[w, i - 1] <= -1 && delta.Wei[w, i] == 1))
      D.Wei[w, i] <- 2 * sum(delta.Wei[w, 1:i]) - i
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
  for (t in 1:length(t.atkin)) {
    SB.Atkin[h, , t] <- 2 * cumsum(G.Atkin[t, ]) / 1:n - 1
    L.Atkin[h, , t] <- D.Atkin[t, ]^2 / 1:n
  }
  for (w in 1:length(wei.t)) {
    SB.Wei[h, , w] <- 2 * cumsum(G.Wei[w, ]) / 1:n - 1
    L.Wei[h, , w] <- D.Wei[w, ]^2 / 1:n
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
Loss.Atkinson <- apply(L.Atkin, c(2, 3), mean)
Sel.Bias.Atkinson <- apply(SB.Atkin, c(2, 3), mean)
Loss.Wei <- apply(L.Wei, c(2, 3), mean)
Sel.Bias.Wei <- apply(SB.Wei, c(2, 3), mean)
Loss.Dom <- colMeans(L.Dom[, , 1])
Sel.Bias.Dom <- colMeans(SB.Dom[, , 1])

# Labels and plot
lab.s <- c(
  paste("Efron(", round(p, 2), ")", sep = ""),
  paste("ABCD(", round(a, 2), ")", sep = ""),
  paste("Atkinson(", t.atkin, ")", sep = ""),
  paste("Wei(", wei.t, ")", sep = ""),
  "Dominant"
)

dfL <- data.frame(
  x = rep(1:n, length(lab.s)),
  y = c(as.vector(Loss.Efron), as.vector(Loss.ABCD), as.vector(Loss.Atkinson), as.vector(Loss.Wei), Loss.Dom),
  des = factor(rep(lab.s, each = n))
)

dfB <- data.frame(
  x = rep(1:n, length(lab.s)),
  y = c(as.vector(Sel.Bias.Efron), as.vector(Sel.Bias.ABCD), as.vector(Sel.Bias.Atkinson), as.vector(Sel.Bias.Wei), Sel.Bias.Dom),
  des = factor(rep(lab.s, each = n))
)

plotL <- ggplot(dfL, aes(x = x, y = y)) +
  geom_line(aes(linetype = des, color = des), size = 1) +
  labs(x = "n", y = "Loss") +
  theme(legend.position = 'bottom') +
  coord_cartesian(xlim = c(15, n))

plotB <- ggplot(dfB, aes(x = x, y = y)) +
  geom_line(aes(linetype = des, color = des), size = 1) +
  labs(x = "n", y = "Selection Bias") +
  theme(legend.position = 'bottom') +
  coord_cartesian(xlim = c(15, n))

# Show plots
windows()
ggarrange(plotL, plotB, nrow = 2, common.legend = TRUE, legend = "bottom")

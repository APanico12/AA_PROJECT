#### Allocation Adaptive Designs ####
rm(list = ls())
# installing/loading the packages:
if(!require(ggplot2)) {
  install.packages("ggplot2"); require(ggplot2)}
if(!require(ggpubr)) {
  install.packages("ggpubr"); require(ggpubr)}
if(!require(progress)) {
  install.packages("progress"); require(progress)}
ptm = proc.time()

# Settings ####
n <- 200 # n. patients
M <- 5000 # MC replicates 
p <- 2/3 # Efron
a <- c(3, 1, 1/3) # ABCD
t.atkin <- 1:2 # Atkinson
set.seed(1)

# create a progress bar
pb <- progress_bar$new(
  format = " [:bar] :percent eta: :eta",
  total = M, clear = FALSE, width= 60)

# empty boxes for delta, G, D, SB and Loss
delta.Efron <- G.Efron <- D.Efron <- matrix(NA, ncol = n, nrow = length(p))
delta.ABCD <- G.ABCD <- D.ABCD <- matrix(NA, ncol = n, nrow = length(a))
delta.Atkin <- G.Atkin <- D.Atkin <- matrix(NA, ncol = n, nrow = length(t.atkin))
SB.Efron <- L.Efron <- array(NA, dim = c(M, n, length(p)))
SB.ABCD <- L.ABCD <- array(NA, dim = c(M, n, length(a)))
SB.Atkin <- L.Atkin <- array(NA, dim = c(M, n, length(t.atkin)))

# MC replications ####
for (h in 1:M) {
  # Initialization: 1st patient  
  delta.Efron[, 1] <- rbinom(length(p), 1, .5)
  delta.ABCD[, 1] <- rbinom(length(a), 1, .5)
  delta.Atkin[, 1] <- rbinom(length(t.atkin), 1, .5)
  G.Efron[, 1] <- rbinom(length(p), 1, .5)
  G.ABCD[, 1] <- rbinom(length(a), 1, .5)
  G.Atkin[, 1] <- rbinom(length(t.atkin), 1, .5)
  D.Efron[, 1] <- 2 * delta.Efron[, 1] - 1
  D.ABCD[, 1] <-  2 * delta.ABCD[, 1] - 1
  D.Atkin[, 1] <- 2 * delta.Atkin[, 1] - 1
  SB.Efron[h , 1, ] <- 2 * G.Efron[, 1] - 1
  SB.ABCD[h , 1, ] <- 2 * G.ABCD[, 1] - 1
  SB.Atkin[h, 1, ] <- 2 * G.Atkin[, 1] - 1
  L.Efron[h , 1, ] <- D.Efron[, 1] ^ 2
  L.ABCD[h , 1, ] <- D.ABCD[, 1] ^ 2
  L.Atkin[h, 1, ] <- D.Atkin[, 1] ^ 2
  
  # start the trial
  for (i in 2:n) {
    # Efron (slide 12) ####
    for (e in 1:length(p)) {# loop for different values of p
      if (D.Efron[e, i - 1] == 0) { # perfectly balanced: D = 0
        delta.Efron[e, i] <- rbinom(1, 1, .5) # p = 0.5
        G.Efron[e, i] <- rbinom(1, 1, .5)
      } else if (D.Efron[e, i - 1] >= 1) { # More A's than B's: D >= 1
        delta.Efron[e, i] <- rbinom(1, 1, 1 - p[e]) # assign A with probability = 1 - p
        if (delta.Efron[e, i] == 0) { # predictability
          G.Efron[e, i] <- 1 # pick the under-represented treatment (correct guess)
        } else {
          G.Efron[e, i] <- 0
        }
      } else if (D.Efron[e, i - 1] <= -1) { # More B's than A's: D <= 1
        delta.Efron[e, i] <- rbinom(1, 1, p[e]) # assign A with probability = p
        if (delta.Efron[e , i] == 1) { # predictability
          G.Efron[e, i] <- 1 # pick the under-represented treatment (correct guess)
        } else {
          G.Efron[e, i] <- 0
        }
      }
    }
    
    # Atkinson: t = 1,2 (slide 13) ####
    for (t in 1:length(t.atkin)) {
      if (D.Atkin[t, i - 1] == 0) {# perfectly balanced: D = 0
        delta.Atkin[t, i] <- rbinom(1, 1, .5) # assign A with probability 0.5
        G.Atkin[t, i] <- rbinom(1, 1, .5)
      } else if (D.Atkin[t, i - 1] >= 1) { # More A's than B's: D >= 1
        delta.Atkin[t, i] <-
          rbinom(1, 1, (sum(1 - delta.Atkin[t, 1:(i - 1)]) ^ t.atkin[t])  /
                   (sum(delta.Atkin[t, 1:(i - 1)]) ^ t.atkin[t] + sum(1 - delta.Atkin[t, 1:(i - 1)]) ^ t.atkin[t]))
        if (delta.Atkin[t, i] == 0) { # predictability
          G.Atkin[t, i] <- 1 # pick the under-represented treatment (correct guess)
        } else {
          G.Atkin[t, i] <- 0
        }
      } else if (D.Atkin[t, i - 1] <= -1) { # More B's than A's: D <= 1
        delta.Atkin[t, i] <-
          rbinom(1, 1, (sum(1 - delta.Atkin[t, 1:(i - 1)]) ^ t.atkin[t])  /
                   (sum(delta.Atkin[t, 1:(i - 1)]) ^ t.atkin[t] + sum(1 - delta.Atkin[t, 1:(i - 1)]) ^ t.atkin[t]))
        if (delta.Atkin[t, i] == 1) {
          # predictability
          G.Atkin[t, i] <- 1 # pick the under-represented treatment (correct guess)
        } else {
          G.Atkin[t, i] <- 0
        }
      }
    }
    
    # ABCD (slide 14) ####
    for (j in 1:length(a)) {# loop for different values of a
      if (D.ABCD[j, i - 1] == 0) { # perfectly balanced: D = 0
        delta.ABCD[j, i] <- rbinom(1, 1, .5) # assign A with probability 0.5
        G.ABCD[j, i] <- rbinom(1, 1, .5)
      } else if (D.ABCD[j, i - 1] >= 1) { # More A's than B's: D >= 1
        delta.ABCD[j, i] <- rbinom(1, 1, 1 / (D.ABCD[j, i - 1] ^ a[j] + 1))
        if (delta.ABCD[j, i] == 0) { # predictability
          G.ABCD[j, i] <- 1 # pick the under-represented treatment (correct guess)
        } else {
          G.ABCD[j, i] <- 0
        }
      } else if (D.ABCD[j, i - 1] <= -1) { # More B's than A's: D <= 1
        delta.ABCD[j, i] <-
          rbinom(1, 1, abs(D.ABCD[j, i - 1]) ^ a[j] / (abs(D.ABCD[j, i - 1]) ^ a[j] + 1))
        if (delta.ABCD[j, i] == 1) { # predictability
          G.ABCD[j, i] <- 1 # pick the under-represented treatment (correct guess)
        } else {
          G.ABCD[j, i] <- 0
        }
      }
    }
    
    
    # Calculate Imbalance D_n
    for (e in 1:length(p)) {
      D.Efron[e, i] <- 2 * sum(delta.Efron[e, 1:i]) - i
    }
    
    for (j in 1:length(a)) {
      D.ABCD[j, i] <- 2 * sum(delta.ABCD[j, 1:i]) - i
    }
    
    for (t in 1:length(t.atkin)) {
      D.Atkin[t, i] <- 2 * sum(delta.Atkin[t, 1:i]) - i
    }
    
  } # end trial
  
  # Calculate Loss and Selection Bias
  for (e in 1:length(p)) {
    SB.Efron[h, , e] <-
      2 * cumsum(G.Efron[e,]) / 1:n - 1
    L.Efron[h, , e] <-
      D.Efron[e,] ^ 2 / 1:n
  }
  
  for (j in 1:length(a)) {
    SB.ABCD[h, , j] <-
      2 * cumsum(G.ABCD[j,]) / 1:n - 1
    L.ABCD[h, , j] <-
      D.ABCD[j,] ^ 2 / 1:n
  }
  for (t in 1:length(t.atkin)) {
    SB.Atkin[h, , t] <- 2 * cumsum(G.Atkin[t,]) / 1:n - 1
    L.Atkin[h, , t] <- D.Atkin[t,] ^ 2 / 1:n
  }
  
  pb$tick()
  #print(h)
}# end MC
elapsed.time = proc.time() - ptm
elapsed.time

# Averaging MC results ####
Loss.Efron <- Sel.Bias.Efron <- matrix(NA, nrow = n, ncol = length(p))
Loss.ABCD <- Sel.Bias.ABCD <- matrix(NA, nrow = n, ncol = length(a))
Loss.Atkinson <- Sel.Bias.Atkinson <- matrix(NA, nrow = n, ncol = length(t.atkin))

for (e in 1:length(p)) {
  Loss.Efron[, e] <-
    apply(L.Efron[, , e], 2, mean)
  Sel.Bias.Efron[, e] <-
    apply(SB.Efron[, , e], 2, mean)
}

for (j in 1:length(a)) {
  Loss.ABCD[, j] <-
    apply(L.ABCD[, , j], 2, mean)
  Sel.Bias.ABCD[, j] <-
    apply(SB.ABCD[, , j], 2, mean)
}
for (t in 1:length(t.atkin)) {
  Loss.Atkinson[, t] <- apply(L.Atkin[, , t], 2, mean)
  Sel.Bias.Atkinson[, t] <- apply(SB.Atkin[, ,t], 2, mean)
}

# Plots ####

#plot Loss
lab.s <- c(paste("Efron(", round(p,2), ")", sep = ""), paste("ABCD(", round(a,2), ")", sep = ""), paste("Atkinson(", t.atkin, ")", sep = ""))
dfL = data.frame(x = rep(1:n, length(p) + length(a) + length(t.atkin)), y = c(as.vector(Loss.Efron), as.vector(Loss.ABCD), as.vector(Loss.Atkinson)), des = factor(rep(lab.s, each = n)))
plotL = ggplot(dfL, aes(x = x, y = y)) + 
  geom_line(aes(linetype = des, col = des), size = 1) + 
  xlab("n") + ylab("Loss") + 
  theme(legend.position = 'bottom', legend.key.width = unit(0.9, 'cm')) +
  scale_color_discrete(name = "Design:") + scale_linetype_discrete(name = "Design:") +
  coord_cartesian(xlim = c(15, n))

#plot SB
dfB = data.frame(x = rep(1:n, length(p) + length(a) + length(t.atkin)), y = c(as.vector(Sel.Bias.Efron), as.vector(Sel.Bias.ABCD), as.vector(Sel.Bias.Atkinson)), des = factor(rep(lab.s, each = n)))
plotB = ggplot(dfB, aes(x = x, y = y)) + 
  geom_line(aes(linetype = des, col = des), size = 1) + 
  xlab("n") + ylab("Bias") + 
  theme(legend.position = 'bottom', legend.key.width = unit(0.9, 'cm')) +
  scale_color_discrete(name = "Design:") + scale_linetype_discrete(name = "Design:") +
  coord_cartesian(xlim = c(15, n))

windows()
ggarrange(plotL, plotB, ncol = 1, nrow = 2, common.legend = T, legend = "bottom")

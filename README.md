# README: Allocation Adaptive Design and Biased Coin Design (BCD)

## Overview

This presentation explores **allocation adaptive designs** in clinical trials, focusing on how to allocate patients to two treatment groups, A and B, in a way that balances the goals of:

- Maintaining **randomization** to reduce bias and ensure valid inference.
- Achieving **balance** in the number of patients allocated to each treatment to maximize statistical power and precision.

We examine the conflict between perfect balance and randomness, and introduce the **Biased Coin Design (BCD)** as a practical compromise. The BCD adaptively adjusts allocation probabilities to balance treatments while preserving some randomness.

---

## What We Are Doing

- **Imbalance and its growth:**  
  We study the imbalance measure $$D_n$$, defined as the difference in counts of patients allocated to treatments A and B after $$n$$ assignments:
  $D_n = 2 \sum_{k=1}^n \delta_k - n,$
  where $$\delta_k = 1$$ if the $$k$$-th patient is assigned to treatment A, otherwise 0.  
  Tracking $$D_n$$ helps us understand how well the design balances groups over time.

- **Balance and randomness trade-off:**  
  Designs that force perfect balance reduce randomness and increase predictability, which may introduce bias. Designs with full randomization maintain unpredictability but can yield large imbalance. BCDs offer a middle ground.

---

## Measures of Balance: Loss and Selection Bias (SB)

Two key measures help assess the quality of the allocation design:

### 1. Loss $$L_n$$

Loss quantifies the variance inflation due to imbalance:

$$
L_n = \left(\frac{D_n}{\sqrt{n}}\right)^2 = 4n \left(\pi_n - \frac{1}{2}\right)^2,
$$

where $$\pi_n$$ is the proportion of subjects allocated to treatment A after $$n$$ patients.

- **Interpretation:**  
  Loss measures how far the allocation proportion deviates from perfect balance (0.5). Smaller values mean better balance, resulting in more precise treatment effect estimates.

- **Expected loss:**  
  The expected loss is related to the variance of $$D_n$$ as
 
  $\mathbb{E}[L_n] = \frac{\mathrm{Var}(D_n)}{n}.$

### 2. Selection Bias (SB)

Standardized bias measures the predictability of the allocation sequence:

$\mathrm{SB} = \frac{|\mathbb{E}[D_n]|}{\sqrt{\mathrm{Var}(D_n)}}$
- **Interpretation:**  
  A high SB indicates predictable allocations (less randomness), while a low SB implies more random and less predictable assignments.

---

## Aim of Biased Coin Design (BCD)

The BCD is designed to:

- **Adaptively balance treatments:** It uses the current imbalance $$D_n$$ to bias the probability of assigning the next patient, giving a higher chance to the underrepresented treatment.
  
- **Balance randomness and control:**  
  By tuning parameters, BCD can range from pure randomization (max randomness, less balance) to deterministic balancing (max balance, no randomness).

- **Maintain asymptotic balance:**  
  For large sample sizes, BCD allows allocation probabilities to approach 0.5, preserving randomness and reducing predictability.

- **Provide flexibility:**  
  Variants such as Adaptive BCD, Bayesian BCD, and Adjustable/Dominant BCD control the rate at which randomness grows with sample size, addressing limitations of earlier designs.

---

## Where to Find the Code / Implementation

- Sensitivity analysis on $$D_n$$ and $$n$$ can be found in [`simulation_imbalance.R`](simulation_imbalance.R)  
- Simulation about Selection Bias and Loss in [`Simulation_complete.R`](Simulation_complete.R)

---

## Summary

This project/presentation covers:

- The theory behind allocation adaptive designs.
- Quantitative measures of imbalance and randomness.
- The trade-off between balance and randomness.
- Analysis of several **Biased Coin Designs (BCD)** including:
  - The classic Efron's BCD,
  - Adaptive BCD,
  - Bayesian BCD,
  - Adjustable and Dominant BCD variants,
  
  which adaptively control imbalance growth and randomness to maintain trial integrity and power.

- Practical approaches to control imbalance growth and maintain trial integrity.

---

If you have any questions or need further explanations, feel free to reach out!

---

*Prepared for [Modern Design of Experiment]*  
*Date: [22/07/2025]*  
*Authors*

 *Antonio Panico*
   PhD candidate in Industrial Engineering, University of Parma, Italy.
 *Umberto Esposito*
   PhD candidate in Statistical Science, University of Bologna, Italy.
---

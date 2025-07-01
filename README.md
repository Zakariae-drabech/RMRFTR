# RMRFTR: Recursive Markov Random Fields for Change Point Detection

This repository contains the implementation of the **Recursive Markov Random Fields Trend Restore (RMRFTR)** method for detecting multiple change points in noisy piecewise linear (PWL) signals with  discontinuities.

> Zakariae Drabech, 2025.

---

## Overview

Change point detection (CPD) plays a central role in many scientific fields such as climatology, genomics, economics, and signal processing. This method targets signals with a **time-varying mean**, modeled as a **piecewise linear (PWL)** function that includes **abrupt jumps** between segments. The primary objective is to detect both the **number and positions of CPs** from noisy measurements.

The noisy observation model is:

\[
y_j = f_j + b_j, \quad b_j \sim \mathcal{N}(0, \sigma^2)
\]

where \( f_j \) is the true signal, and \( b_j \) is i.i.d. Gaussian noise.

---

## Signal Model

The underlying signal \( f(t) \) is modeled as:

\[
f(t) = \sum_{j=1}^{K+1} \left( \alpha_j + \beta_j t \right) \cdot \mathds{1}_{[\tau_{j-1} + 1, \tau_j]}(t)
\]

- \( \tau_1, \dots, \tau_K \) are the unknown change points  
- \( \alpha_j \), \( \beta_j \) are the intercept and slope of the \( j \)-th segment  
- The signal is **discontinuous** at change points:
  
\[
f(\tau_j^-) + \beta_j \neq f(\tau_j^+), \quad j = 1, \dots, K
\]

This model allows for **slope changes** and **abrupt jumps** between segments.

---

## Key Features

- Detects change points in **discontinuous piecewise linear** signals.
- Assumes **Gaussian noise** with unknown variance.
- Based on a **Markov Random Field (MRF)** formulation.
- Uses a **recursive algorithm** to segment data efficiently.
  
---

## Dependencies

- Python ≥ 3.6  
- `numpy`, `scipy`

# RMRFTR: Recursive Markov Random Fields Trend Restore 

This repository contains the implementation of the **Recursive Markov Random Fields Trend Restore (RMRFTR)** method for detecting multiple change points in noisy piecewise linear (PWL) signals with  discontinuities.

> Zakariae Drabech, 2025.

---

## Overview

Change point detection (CPD) plays a central role in many scientific fields such as climatology, genomics, economics, and signal processing. This method targets signals with a **time-varying mean**, modeled as a **piecewise linear (PWL)** function that includes **abrupt jumps** between segments. The primary objective is to detect both the **number and positions of CPs** from noisy measurements.

The noisy observation model is:

$$y_j = f_j + b_j, \quad b_j \sim \mathcal{N}(0, \sigma^2)$$

where $f_j$ is the true signal, and $b_j$ is i.i.d. Gaussian noise.

---

## Signal Model

The underlying signal $f(t)$ is modeled as:

$$f(t) = \sum_{j=1}^{K+1} ( \alpha_j + \beta_j t ) 1_{[\tau_{j-1} + 1, \tau_j]}(t),$$

- $\tau_1, \dots, \tau_K$ are the unknown change points  
- $\alpha_j$, $\beta_j$ are the intercept and slope of the $j$-th segment  
- The signal is **discontinuous** at change points:
  
$$f(\tau_j)+\beta_{j} \neq f(\tau_j+1),\; j=1, \cdots, K.$$

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

## Quick Start

from RecMRFTR import RMRFTR
import numpy as np
import matplotlib.pyplot as plt

# Step 1: Generate a synthetic PWL signal with jumps
def f(x):
    if x <= 100:
        return x / 50                          
    elif x <= 200:
        return 4 + (x - 100) / 20              
    else:
        return 6 - (x - 200) / 100             

        
n = 300
t=np.arange(n)          
ff=np.vectorize(f)
signal=ff(t)

# Step 2: Add Gaussian noise
np.random.seed(0)
sigma = 1
y = signal + np.random.normal(0, sigma, size=n)

# Step 3: Run RMRFTR
model = RMRFTR(mu=50, h0=1.5)
model.run(y)
restored = model.x_mrf

# Step 4: Plot results
plt.figure(figsize=(12, 4))
plt.plot(y, 'ko',label="Noisy Signal", alpha=0.6)
plt.plot(signal,'g' , label="True Signal",  linewidth=3)
plt.plot(restored,'r', label="Restored Signal", linewidth=3)
for cp in model.ChangePoints:
    plt.axvline(cp, color="red", linestyle=":", label="Detected CP" if cp == model.ChangePoints[0] else "")
plt.legend()
plt.title("RMRFTR: Change Point Detection in Noisy PWL Signal")
plt.tight_layout()
plt.show()

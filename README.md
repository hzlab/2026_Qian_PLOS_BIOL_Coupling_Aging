# Structure–Function Coupling/Decoupling (Graph-Spectral Framework)

This repository contains code used in:

**Xing Qian et al. (2026)**  
Altered salience network structure–function integration underlies the decline in cognitive flexibility during aging  
PLOS Biology

---

## Overview

This implementation is based on:

Medaglia et al. (2018), Nature Human Behaviour

### Key difference

Instead of selecting a fixed number of eigenmodes, we use a **variance-based criterion** to determine:

- Coupling (alignment; smooth components)
- Decoupling (liberality; rough components)

This approach:
- Is data-driven
- Is robust across parcellations
- Avoids arbitrary selection of k

---

## Methods

### 1. Component selection (choose_k.m)

Steps:

1. Divide eigenvalues into:
   - Smooth (positive eigenvalues)
   - Rough (negative eigenvalues)

2. For each group:
   - Compute total energy = sum of absolute eigenvalues
   - Sort eigenvalues (descending)
   - Compute cumulative energy

3. Define kA and kL as:
   - Minimum number of eigenvalues explaining 95% energy

4. Average across time → KA, KL

---

### 2. Signal decomposition (signal_decomposition.m)

Structural matrix:

A = V Λ V^T

- V: eigenvectors
- Λ: eigenvalues

Interpretation:
- Smooth modes → coherent activity
- Rough modes → abrupt variation

---

### Graph Fourier Transform

x_hat = V^T x

Reconstruction:

x = V x_hat

---

### Coupling / Decoupling

- Coupling = projection onto smooth modes
- Decoupling = projection onto rough modes

Use KA and KL to select components.

---

## Notes

- This is NOT SC–FC correlation
- It is graph spectral decomposition

---

## Citation

Please cite:
- Qian et al., 2026 (PLOS Biology)
- Medaglia et al., 2018

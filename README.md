# Exploring the effect of random oscillator heterogeneity on network synchronization

This repository is part of the paper: Y. Zhang, J. L. Ocampo-Espindola, I. Z. Kiss and A. E. Motter, _Random heterogeneity outperforms design in network synchronization_.

We hope the scripts and data included here will facilitate others in reproducing our results.
The repository includes data analysis protocol for time series measured from the experiments, code for simulating delay-coupled Stuart-Landau oscillators in the presence of oscillator heterogeneity, and code for analyzing the stability of synchronization states.

All scripts below can be run independent of each other.

1. `analysis.py` and `analysis.ipynb`

  _This script analyzes how the level of synchrony changes with oscillator heterogeneity in our electrochemical experiments._

2. `analysis2.py` and `analysis2.ipynb`

  _This script compares homogeneous and heterogeneous oscillators in terms of the measured oscillator heterogeneity (when uncoupled) and the time-averaged synchronization error (when coupled)._

3. `plot_trj.py` and `plot_trj.ipynb`

  _Plotting the time series of oscillators measured from electro-chemical experiments._

4. `stuart_landau_delay_dynamics.m`

  _Simulating delay-coupled Stuart-Landau oscillators using MATLAB dde23._

5. `stability_code.m`

  _MATLAB code calculating the maximal transverse Lyapunov exponent of the synchronization state by solving (transcedental) characteristic equations._

6. data

  _This folder contains the time series measured from the experiments._

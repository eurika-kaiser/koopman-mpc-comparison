# Model predictive control of nonlinear systems using Koopman models

This repository contains Matlab code used for generating figures in Fig. 6.6 of the paper "Modern Koopman Theory for Dynamical Systems" by S. L. Brunton, M. Budisic, E. Kaiser, and J. N. Kutz. It's preprint is available on arXiv under: https://arxiv.org/abs/2102.12086

The code compares Koopman-based MPC using DMDc, eDMDc, and time-delay models for several dynamical systems. The code is an adaption/extension of the code associated with the paper "Linear predictors for nonlinear dynamical systems: Koopman operator meets model predictive control", Automatica 2018, by Milan Korda and Igor Mezic (available under: https://github.com/MilanKorda/KoopmanMPC).

### Installation

This package requires the installation of qpOASES and source code of Koopman-MPC (both available under: https://github.com/MilanKorda/KoopmanMPC).

Follow the installation recommendations for qpOASES and add the KoopmanMPC/Resources folder to the Matlab search path, e.g. by executing addpath('./KoopmanMPC/Resources').

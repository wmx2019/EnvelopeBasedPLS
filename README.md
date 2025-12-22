# Envelope-based Partial Least Squares in Functional Regression

This repository contains R code for the paper:

**"Envelope-based Partial Least Squares in Functional Regression"**  
Minxuan Wu, Joseph Antonelli, Zhihua Su  
Accepted pending minor revision at *Journal of Multivariate Analysis*  
[arXiv:2505.14876](https://arxiv.org/abs/2505.14876)

## Overview

This paper addresses the gap in functional PLS literature by providing the first methods with provable theoretical guarantees within the functional predictor envelope framework. We develop two methods:
- **FEPLS**: Functional Envelope-based PLS for scalar and functional responses
- **GFEPLS**: Generalized Functional Envelope-based PLS for binary and other response distributions

Both methods achieve root-n consistency for finite-rank coefficient functions, with FEPLS remaining consistent even when rank grows with sample size. We establish asymptotic normality and provide pointwise confidence and prediction intervals.

## Repository Structure

```
JMA/
├── FEPLS.R                      # Main implementation of FEPLS and GFEPLS methods
├── Basic_Functions.R            # Helper functions
├── DA_binary.R                  # Data analysis: binary response
├── DA_functional.R              # Data analysis: functional response
├── DS_binary.R                  # Data simulation: binary response
├── DS_scalar.R                  # Data simulation: scalar response
├── DS_functional.R              # Data simulation: functional response
└── Functional_Response_Data/    # Real data for functional response analysis
```

## Files Description

### Core Implementation
- **FEPLS.R**: Contains the main implementation of both FEPLS (for scalar/functional responses) and GFEPLS (for binary responses)
- **Basic_Functions.R**: Supporting functions used across different methods

### Data Analysis (DA)
Scripts for analyzing real data with different response types:
- **DA_binary.R**: Binary response with functional predictor
- **DA_functional.R**: Functional response with functional predictor

### Data Simulation (DS)
Scripts for simulation studies with different response types:
- **DS_binary.R**: Binary response simulation
- **DS_scalar.R**: Scalar response simulation  
- **DS_functional.R**: Functional response simulation

## Usage

To reproduce the results from the paper:

**Run simulations:**
```r
# For scalar response
source("JMA/DS_scalar.R")

# For binary response
source("JMA/DS_binary.R")

# For functional response
source("JMA/DS_functional.R")
```

**Run real data analysis:**
```r
# For binary response analysis
source("JMA/DA_binary.R")

# For functional response analysis
source("JMA/DA_functional.R")
```
# QFCV

Quantile-based forward cross-validation for predictive intervals for prediction error in time series forecasting. 

This method was introduced in the paper Uncertainty intervals for Prediction Errors in Time Series
Forecasting by Hui Xu, Song Mei, Stephen Bates, Jonathan Taylor, and Robert Tibshirani. 

The repository contains four folders, corresponding code to reproduce plots in different sections of the paper. 

### Motivation 
- **helper.R**: contains helper functions including data generation, model fitting, FCV, and QFCV
- **compare.R**: contains code to compare stochastic and expected test errors, along with naive FCV and QFCV intervals for 50 simulation instances. This reproduces Figure 1 in the paper.
- **histogram.R**: contains code to reproduce Figure 10 in the paper.
- **data**
  
### QFCV simulation

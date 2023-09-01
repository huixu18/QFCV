# QFCV

Quantile-based forward cross-validation for predictive intervals for prediction error in time series forecasting. 

This method was introduced in the paper Uncertainty intervals for Prediction Errors in Time Series
Forecasting by Hui Xu, Song Mei, Stephen Bates, Jonathan Taylor, and Robert Tibshirani. 

The repository contains four folders, corresponding code to reproduce plots in different sections of the paper listed as follows. 

### Motivation 
- **helper.R**: contains helper functions including data generation, model fitting, FCV, and QFCV
- **compare.R**: contains code to compare stochastic and expected test errors, along with naive FCV and QFCV intervals for 50 simulation instances. This reproduces Figure 1 in the paper.
- **histogram.R**: contains code to reproduce Figure 10 in the paper.
  
### QFCV simulation
- **helper.R**: contains helper functions including data generation, model fitting, FCV, and QFCV (with different variants)
- **QFCV_variants_5_AR1.R**: contains code to numerically compare QFCV variants for test size 5 and ARMA(1,0) noise (Figure 6(a))
- **QFCV_variants_5_AR1_MA20.R**: contains code to numerically compare QFCV variants for test size 5 and ARMA(1,20) noise (Figure 6(b))
- **QFCV_variants_20_AR1.R**: contains code to numerically compare QFCV variants for test size 20 and ARMA(1,0) noise (Figure 6(c))
- **QFCV_variants_20_AR1_MA20.R**: contains code to numerically compare QFCV variants for test size 20 and ARMA(1,20) noise (Figure 6(d))

### FCV
- **helper.R**: contains helper functions including data generation, model fitting, FCV with its variants, and QFCV(1)
- **FCV_variants_5_AR1.R**: contains code to numerically compare FCV variants with QFCV(1) for test size 5 and ARMA(1,0) noise (Figure 11(a))
- **FCV_variants_5_AR1_MA20.R**: contains code to numerically compare FCV variants with QFCV(1) for test size 5 and ARMA(1,20) noise (Figure 11(b))
- **FCV_variants_20_AR1.R**: contains code to numerically compare FCV variants with QFCV(1) for test size 20 and ARMA(1,0) noise (Figure 11(c))
- **FCV_variants_20_AR1_MA20.R**: contains code to numerically compare FCV variants with QFCV(1) for test size 20 and ARMA(1,20) noise (Figure 11(d))

### AQFCV
- **helper.R**: contains helper functions including data generation, model fitting, FCV, and AQFCV
- **coverage.R**: contains code to generate 

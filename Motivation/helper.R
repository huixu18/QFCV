##############################################
# Helper functions for generating motivation plots
##############################################

library(glmnet)
library(quantreg)
require(doMC)
registerDoMC(cores = 4)

###########################
# Data generation 
###########################
# Params: n = number of training samples, p = number of features, n_test = number of test samples, alpha = autocorrelation
# Output: x = training data matrix, y = training outcome vector, x_test = test data matrix, y_test = test outcome vector 
data_gen <- function(n, p, n_test, alpha = 0){
  beta = rep(0, p)
  beta[1:4] = 1
  x = matrix(rnorm(n * p), nrow = n, ncol = p)
  x_test = matrix(rnorm(n_test * p), nrow = n_test, ncol = p)
  noise = arima.sim(list(order = c(1,0,0), ar = alpha), n = n + n_test)
  y = x %*% beta + noise[1:n]
  y_test = x_test %*% beta + noise[(n+1):(n+n_test)]
  return(list(x = x, y = y, x_test = x_test, y_test = y_test))
}

###########################
# Model fitting and evaluation using glmnet
###########################
# Params: x_train = training data matrix, y_train = training outcome vector, x_test = test data matrix, y_test = test outcome vector
# Output: err_list = vector of squared loss on test data
fit <- function(x_train, y_train, x_test, y_test){
  a=glmnet(x_train, y_train, family="gaussian", parallel = TRUE, lambda = 0.05)
  y_hat=predict(a,x_test)
  err_list = (y_hat - y_test)^2
  return(err_list)
}

###########################
# Forward cross-validation (FCV)
###########################
# Params: n = number of samples, x = data matrix, y = outcome vector, w = window size, w_train = window size for training, lag = FCV lag 
# Output: err_val = vector of validation errors obtained from FCV 
fcv <- function(n, x, y, w, w_train, lag){
  start = 1
  end = w 
  nfold = floor((n-w)/lag) + 1
  err_val = c()
  for(fold in 1:nfold){
    split = start + w_train - 1
    err_val = c(err_val, mean(fit(x[start:split,], y[start:split], x[(split+1):end,], y[(split+1):end])))
    start = start + lag
    end = end + lag 
  }
  return(err_val)
}

###########################
# Quantile-based forward cross-validation (QFCV)
###########################
# Params: n = number of samples, x = data matrix, y = outcome vector, w = window size, w_train = window size for training, w_val = window size for validation, lag = FCV lag
# Output: err_val_list = vector of validation errors, err_test_list = vector of test errors, err_hat_lo = lower end of predictive interval, err_hat_hi = upper end of predictive interval, err_hat = point estimator for prediction error 
fcv_cal <- function(n, x, y, w, w_train, w_val, lag){
  w_test = w - w_train - w_val 
  start = 1
  end = w 
  nfold = floor((n-w)/lag) + 1
  err_val = c()
  err_test = c()
  for(fold in 1:nfold){
    train_end = start + w_train - 1
    val_start = start + w_train 
    val_end = start + w_train + w_val - 1
    test_start = start + w_train + w_val
    err_val = c(err_val, mean(fit(x[start:train_end,], y[start:train_end], x[val_start:val_end,], y[val_start:val_end])))
    err_test = c(err_test, mean(fit(x[(val_end - w_train + 1):val_end,], y[(val_end - w_train + 1):val_end], x[test_start:end,], y[test_start:end])))
    start = start + lag
    end = end + lag 
  }

  val = fit(x[(n-w_train-w_val+1):(n-w_val),], y[(n-w_train-w_val+1):(n-w_val)], x[(n-w_val+1):n,], y[(n-w_val+1):n])
  point_val = mean(val)
  y_reg = err_test
  x_reg = err_val
  lr = lm(y_reg ~ x_reg)
  err_hat = predict(lr, data.frame(x_reg = point_val))
  
  rqdata <- data.frame(err_test = err_test, err_val = err_val)
  rqfit_lo <- rq(err_test ~ err_val, tau = 0.05, data = rqdata)
  err_hat_lo = predict(rqfit_lo, newdata = data.frame(err_val = c(point_val)))
  rqfit_hi <- rq(err_test ~ err_val, tau = 0.95, data = rqdata)
  err_hat_hi = predict(rqfit_hi, newdata = data.frame(err_val = c(point_val)))
  
  return(list(err_val_list = err_val, err_test_list = err_test,
              err_hat_lo = err_hat_lo, err_hat_hi = err_hat_hi, err_hat = err_hat))
}


##############################################
# Helper functions for generating motivation plots
##############################################

library(glmnet)
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
# Params: n = number of training samples, x = data matrix, y = outcome vector, w = window size, w_train = window size for training, lag = FCV lag 
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



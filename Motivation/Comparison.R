library(glmnet)
library(quantreg)
library(ggplot2)
source("utility.R")

fcv_cal <- function(n, x, y, w, w_train, w_val, lag){
  w_test = w - w_train - w_val 
  start = 1
  end = w 
  nfold = floor((n-w)/lag) + 1
  err_val = c()
  err_test = c()
  err_val_mid = c()
  err_val_matrix = matrix(NA, nrow = nfold, ncol = w_val)
  for(fold in 1:nfold){
    train_end = start + w_train - 1
    val_start = start + w_train 
    val_end = start + w_train + w_val - 1
    test_start = start + w_train + w_val
    err_val = c(err_val, mean(fit(x[start:train_end,], y[start:train_end], x[val_start:val_end,], y[val_start:val_end])))
    val_end_mid = start + w_train + floor(w_val/2) - 1
    err_val_mid = c(err_val_mid, mean(fit(x[start:train_end,], y[start:train_end], x[val_start:val_end_mid,], y[val_start:val_end_mid])))
    err_val_matrix[fold,] = fit(x[start:train_end,], y[start:train_end], x[val_start:val_end,], y[val_start:val_end])
    err_test = c(err_test, mean(fit(x[(val_end - w_train + 1):val_end,], y[(val_end - w_train + 1):val_end], x[test_start:end,], y[test_start:end])))
    # err_test = c(err_test, mean(fit(x[start:val_end,], y[start:val_end], x[test_start:end,], y[test_start:end])))
    start = start + lag
    end = end + lag 
  }
  ## calibration without covariates 
  err_hat_0 = mean(err_test)
  err_hat_0_lo = quantile(err_test, 0.05)
  err_hat_0_hi = quantile(err_test, 0.95)
  
  ## calibration using point validation error
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
  
  ## calibration using two validation errors 
  lr_mid = lm(err_test ~ err_val + err_val_mid)
  err_hat_mid = predict(lr_mid, data.frame(err_val = point_val, err_val_mid = mean(val[1:floor(w_val/2)])))
  rqdata_mid <- data.frame(err_test = err_test, err_val = err_val, err_val_mid = err_val_mid)
  rqfit_mid_lo <- rq(err_test ~., tau = 0.05, data = rqdata_mid)
  err_hat_mid_lo = predict(rqfit_mid_lo, newdata = data.frame(err_val = c(point_val), err_val_mid = mean(val[1:floor(w_val/2)])))
  rqfit_mid_hi <- rq(err_test ~., tau = 0.95, data = rqdata_mid)
  err_hat_mid_hi = predict(rqfit_mid_hi, newdata = data.frame(err_val = c(point_val), err_val_mid = mean(val[1:floor(w_val/2)])))
  
  ## calibration using list of validation errors 
  rqdata_full = as.data.frame(err_val_matrix)
  rqdata_full = cbind(err_test, rqdata_full)
  lr_full = lm(err_test ~., data = rqdata_full)
  err_hat_full = predict(lr_full, as.data.frame(t(val)))
  rqfit_full_lo <- rq(err_test ~., tau = 0.05, data = rqdata_full)
  err_hat_full_lo <- predict(rqfit_full_lo, newdata = as.data.frame(t(val))) 
  rqfit_full_hi <- rq(err_test ~., tau = 0.95, data = rqdata_full)
  err_hat_full_hi <- predict(rqfit_full_hi, newdata = as.data.frame(t(val))) 
  
  return(list(err_val_list = err_val, err_test_list = err_test, err_hat_0 = err_hat_0, err_hat_0_lo = err_hat_0_lo, err_hat_0_hi = err_hat_0_hi, 
              err_hat_lo = err_hat_lo, err_hat_hi = err_hat_hi, err_hat = err_hat, 
              err_hat_mid = err_hat_mid, err_hat_mid_lo = err_hat_mid_lo, err_hat_mid_hi = err_hat_mid_hi,
              err_hat_full = err_hat_full, err_hat_full_lo = err_hat_full_lo, err_hat_full_hi = err_hat_full_hi))
}

nsim = 500
n = 1000
n_test = 20
p = 20
w = 80
w_train = 40
w_val = 20
lag = 5
cc = TRUE
# alpha_list = c(0.1, 0.3, 0.5, 0.7, 0.9)
alpha_list = c(0.5)
for(alpha in alpha_list){
  err_test = err_fcv_cal = lo_fcv_cal = hi_fcv_cal = err_fcv_cal_mid = lo_fcv_cal_mid = hi_fcv_cal_mid = err_fcv_cal_full = 
    lo_fcv_cal_full = hi_fcv_cal_full = err_fcv_cal_0 = lo_fcv_cal_0 = hi_fcv_cal_0 = err_nfcv = se_nfcv = err_fcv = se_fcv = rep(NA, nsim)
  for(i in 1:nsim){
    data = data_gen(n, p, n_test, alpha = alpha, binary = FALSE)
    cat("Bayes error for iter ", i, " = ", data$err_bayes, fill = TRUE)
    err_test[i] = mean(fit(data$x[(n-w_train+1):n,], data$y[(n-w_train+1):n], data$x_test, data$y_test))
    err_fcv_list = fcv(n, data$x, data$y, w_train+n_test, w_train, 5)
    err_fcv[i] = mean(err_fcv_list)
    if(cc == FALSE){ 
      se_fcv[i] = sd(err_fcv_list)/sqrt(length(err_fcv_list))
    } 
    else{
      var_fcv = var(err_fcv_list)
      m = floor(length(err_fcv_list) / 50) 
      for(s in 1:m){
        var_fcv = var_fcv + 2 * (1 - s/length(err_fcv_list)) * cov(err_fcv_list[1:(length(err_fcv_list)-s)], err_fcv_list[(s+1):(length(err_fcv_list))])
      }
      if(var_fcv < 0){
        var_fcv = var(err_fcv_list)
      }
      se_fcv[i] = sqrt(var_fcv/length(err_fcv_list))
    }
    cat("true test: ", err_test[i], "quantile: ", c(err_fcv_cal[i], hi_fcv_cal[i] - lo_fcv_cal[i]), "se_fcv: ", se_fcv[i], "m: ", m, fill = TRUE)
  }
  
  out = list(err_test = err_test, err_fcv = err_fcv, se_fcv = se_fcv,
             err_fcv_cal = err_fcv_cal, lo_fcv_cal = lo_fcv_cal, hi_fcv_cal = hi_fcv_cal,
             err_fcv_cal_0 = err_fcv_cal_0, lo_fcv_cal_0 = lo_fcv_cal_0, hi_fcv_cal_0 = hi_fcv_cal_0,
             err_fcv_cal_mid = err_fcv_cal_mid, lo_fcv_cal_mid = lo_fcv_cal_mid, hi_fcv_cal_mid = hi_fcv_cal_mid,
             err_fcv_cal_full = err_fcv_cal_full, lo_fcv_cal_full = lo_fcv_cal_full, hi_fcv_cal_full = hi_fcv_cal_full)
  save(out, file = paste0("sparse_linear_20_cc_alpha-", alpha,".RData"))
  print(paste0("Results for alpha-", alpha, " saved to disk."))
}

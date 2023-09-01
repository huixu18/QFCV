##############################################
# Comparison of rolling predictive intervals for ACI-DF wrapped with QFCV and baseline FCV method
##############################################
library(glmnet)
library(quantreg)
library(ggplot2)
library(wesanderson)

source("helper.R")

set.seed(100)
n = 3000
n0 = 500 ## number of historical data points before start training 
n_test = 5
p = 20
w = 50
w_train = 40
w_val = 5
lag = 5
lag_fcv = n_test
gamma = 0.005
gamma_fcv = 10
data = data_gen(n, p, 0, alpha = 0.5, MA = TRUE)
nsim = floor((n-n0)/lag)
err_test = err_qfcv = lo_qfcv = hi_qfcv = lo_qfcv_online = hi_qfcv_online = lo_fcv = hi_fcv = lo_fcv_online = hi_fcv_online = rep(NA, nsim)
fold_start = floor((n0-w)/lag) 
err_val_list = err_test_list = c()
theta = 0
theta_list = c(0)
theta_fcv = 0
theta_fcv_list = c(0)
err_fcv_list_full = fcv(n, data$x, data$y, w_train+n_test, w_train, lag_fcv)
for(i in 1:nsim){
  cat("iter ", i, fill = TRUE)
  split = n0+(i-1) * lag
  err_test[i] = mean(fit(data$x[(split-w_train+1):split,], data$y[(split-w_train+1):split], data$x[(split+1):(split+n_test),], data$y[(split+1):(split+n_test)]))
  fold_fcv = floor((split - w_train - n_test)/lag_fcv) + 1
  err_fcv_list = err_fcv_list_full[1:fold_fcv]
  err_fcv = mean(err_fcv_list)
  se_fcv = sd(err_fcv_list)/sqrt(length(err_fcv_list))
  lo_fcv[i] = err_fcv - qnorm(0.95) * se_fcv 
  hi_fcv[i] = err_fcv + qnorm(0.95) * se_fcv 
  lo_fcv_online[i] = err_fcv - (qnorm(0.95) + theta_fcv) * se_fcv
  hi_fcv_online[i] = err_fcv + (qnorm(0.95) + theta_fcv) * se_fcv
  theta_fcv = theta_fcv + gamma_fcv * (as.numeric(err_test[i] > hi_fcv_online[i] | err_test[i] < lo_fcv_online[i])-0.1)
  result_qfcv = fcv_cal(split, fold_start * (i >1), data$x[1:split,], data$y[1:split], w, w_train, w_val, lag, 0, err_val_list, err_test_list)
  result_qfcv_online = fcv_cal(split, fold_start * (i >1), data$x[1:split,], data$y[1:split], w, w_train, w_val, lag, theta, err_val_list, err_test_list)
  lo_qfcv_online[i] = result_qfcv_online$err_hat_lo
  hi_qfcv_online[i] = result_qfcv_online$err_hat_hi
  
  err_val_list = result_qfcv$err_val_list
  err_test_list = result_qfcv$err_test_list
  
  err_qfcv[i] = result_qfcv$err_hat
  lo_qfcv[i] = result_qfcv$err_hat_lo
  hi_qfcv[i] = result_qfcv$err_hat_hi
  theta = theta + gamma * (as.numeric(err_test[i] > hi_qfcv_online[i] | err_test[i] < lo_qfcv_online[i])-0.1)
  theta_list = c(theta_list, theta)
  theta_fcv_list = c(theta_fcv_list, theta_fcv)
  fold_start = fold_start + 1
  cat("true test: ", err_test[i], "quantile: ", c(lo_qfcv[i], hi_qfcv[i]), fill = TRUE)
}

coverage = as.numeric(err_test < hi_qfcv & err_test > lo_qfcv)
coverage_online = as.numeric(err_test < hi_qfcv_online & err_test > lo_qfcv_online)
coverage_fcv = as.numeric(err_test < hi_fcv & err_test > lo_fcv)
coverage_fcv_online = as.numeric(err_test < hi_fcv_online & err_test > lo_fcv_online)

data_PI_combined = data.frame(time = rep(seq(length(coverage)), 2), method = rep(c("AQFCV", "AFCV"), each = length(coverage)),
                     true = rep(err_test, 2), estimate = rep(err_qfcv, 2), hi = c(hi_qfcv_online, hi_fcv_online), lo = c(lo_qfcv_online, lo_fcv_online),
                     coverage = c(coverage_online, coverage_fcv_online))
ggplot(data_PI_combined) + geom_line(aes(x = time, y = true), size = 0.8) +
  geom_ribbon(aes(x = time, y = estimate, ymin = lo, ymax = hi, col = method, fill = method), alpha = 0.3) + scale_color_manual(values = wes_palette("Chevalier1", n = 2)) + scale_fill_manual(values = wes_palette("Chevalier1", n = 2)) + xlab("Time") + ylab("Test error") + 
  theme_bw() + theme(legend.position="top", text = element_text(family = "Palatino", size = 16))

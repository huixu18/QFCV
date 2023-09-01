##############################################
# Comparison of instance and time average coverage for ACI-DF wrapped with QFCV and baseline FCV method
##############################################
library(ggplot2)
library(ggpubr)
source("helper.R")

###########################
# Run simulation for stationary data
###########################
set.seed(100)
n = 3000
n0 = 500 ## number of historical data points before start training 
n_test = 5
p = 20
w = 50
w_train = 40
w_val = 5
lag = 5
lag_fcv = n_test + w_train
gamma = 0.005
gamma_fcv = 15
rep = 500
nsim = floor((n-n0)/lag) - n_test/lag
cc = FALSE
alpha = 0.5
coverage = coverage_online = coverage_fcv = coverage_fcv_online = matrix(NA, nrow = rep, ncol = nsim)
length = length_online = length_fcv = length_fcv_online = matrix(NA, nrow = rep, ncol = nsim)
for(r in 1:rep){
  cat("iter ", r, fill = TRUE)
  data = data_gen(n, p, 0, alpha = alpha)
  err_test = err_qfcv = lo_qfcv = hi_qfcv = lo_qfcv_online = hi_qfcv_online = lo_fcv = hi_fcv = lo_fcv_online = hi_fcv_online = rep(NA, nsim)
  fold_start = floor((n0-w)/lag) 
  err_val_list = err_test_list = c()
  err_fcv_list = c()
  theta = 0
  theta_list = c(0)
  theta_fcv = 0
  err_fcv_list_full = fcv(n, data$x, data$y, w_train+n_test, w_train, lag_fcv)
  for(i in 1:nsim){
    split = n0+(i-1) * lag
    err_test[i] = mean(fit(data$x[(split-w_train+1):split,], data$y[(split-w_train+1):split], data$x[(split+1):(split+n_test),], data$y[(split+1):(split+n_test)]))
    fold_fcv = floor((split-w_train-n_test)/lag_fcv) + 1
    err_fcv_list = err_fcv_list_full[1:fold_fcv]
    err_fcv = mean(err_fcv_list)
    if(cc == FALSE){ 
      se_fcv = sd(err_fcv_list)/sqrt(length(err_fcv_list))
    } 
    else{
      var_fcv = var(err_fcv_list)
      m = floor(length(err_fcv_list) / 50) 
      for(s in 1:m){
        var_fcv = var_fcv + 2 * (1 - s/length(err_fcv_list)) * cov(err_fcv_list[1:(length(err_fcv_list)-s)], err_fcv_list[(s+1):(length(err_fcv_list))])
      }
      var_fcv = (var_fcv < 0) * var(err_fcv_list) + (var_fcv >= 0) * var_fcv
      se_fcv = sqrt(var_fcv/length(err_fcv_list))
    }
    
    lo_fcv[i] = err_fcv - qnorm(0.95) * se_fcv
    hi_fcv[i] = err_fcv + qnorm(0.95) * se_fcv
    lo_fcv_online[i] = err_fcv - (qnorm(0.95) + theta_fcv) * se_fcv 
    hi_fcv_online[i] = err_fcv + (qnorm(0.95) + theta_fcv) * se_fcv 
    theta_fcv = theta_fcv + gamma_fcv * (as.numeric(err_test[i] > hi_fcv_online[i] | err_test[i] < lo_fcv_online[i])-0.1)
    
    result_qfcv = fcv_cal(split, (fold_start + i - 1) * (i > 1), data$x[1:split,], data$y[1:split], w, w_train, w_val, lag, 0, err_val_list, err_test_list)
    result_qfcv_online = fcv_cal(split, (fold_start + i - 1) * (i > 1), data$x[1:split,], data$y[1:split], w, w_train, w_val, lag, theta, err_val_list, err_test_list)
    lo_qfcv_online[i] = result_qfcv_online$err_hat_lo
    hi_qfcv_online[i] = result_qfcv_online$err_hat_hi
    err_val_list = result_qfcv$err_val_list
    err_test_list = result_qfcv$err_test_list
    err_qfcv[i] = result_qfcv$err_hat
    lo_qfcv[i] = result_qfcv$err_hat_lo
    hi_qfcv[i] = result_qfcv$err_hat_hi
    theta = theta + gamma * (as.numeric(err_test[i] > hi_qfcv_online[i] | err_test[i] < lo_qfcv_online[i])-0.1)
    theta_list = c(theta_list, theta)
  }
  coverage[r,] = as.numeric(err_test < hi_qfcv & err_test > lo_qfcv)
  coverage_online[r,] = as.numeric(err_test < hi_qfcv_online & err_test > lo_qfcv_online)
  coverage_fcv[r,] = as.numeric(err_test < hi_fcv & err_test > lo_fcv)
  coverage_fcv_online[r,] = as.numeric(err_test < hi_fcv_online & err_test > lo_fcv_online)
  length[r,] = hi_qfcv - lo_qfcv
  length_online[r,] = hi_qfcv_online - lo_qfcv_online
  length_fcv[r,] = hi_fcv - lo_fcv
  length_fcv_online[r,] = hi_fcv_online - lo_fcv_online
  cat("true test: ", mean(err_test), "coverage: ", mean(coverage[r,]), "coverage_online: ", mean(coverage_online[r,]), "coverage_fcv: ", mean(coverage_fcv[r,]), "coverage_fcv_online: ", mean(coverage_fcv_online[r,]), fill = TRUE)
}

out = list(coverage = coverage, coverage_online = coverage_online, coverage_fcv = coverage_fcv, coverage_fcv_online = coverage_fcv_online,
           length = length, length_online = length_online, length_fcv = length_fcv, length_fcv_online = length_fcv_online)
save(out, file = paste0("AR1_5_stationary.RData"))

###########################
# Look at results to plot coverage
###########################
load(file = "AR1_5_stationary.RData")
coverage = out$coverage
coverage_online = out$coverage_online
coverage_fcv = out$coverage_fcv
coverage_fcv_online = out$coverage_fcv_online

###########################
# Plot instance-average coverage
###########################
data_coverage = data.frame(coverage = c(apply(coverage, 2, mean), apply(coverage_online, 2, mean)),
                           time = rep(seq(ncol(coverage)), 2), Method = rep(c("QFCV", "AQFCV"), each = ncol(coverage)))
sd = apply(coverage, 2, sd)
p1 <- ggplot(data_coverage) + geom_line(aes(x = time, y = coverage, col = Method)) + theme_bw() + geom_hline(yintercept=0.9, linetype="dashed", size=0.8) + 
  geom_ribbon(aes(x = time, ymin = coverage - qnorm(0.95) * sd/sqrt(500), ymax = coverage + qnorm(0.95) * sd/sqrt(500), fill = Method), alpha = 0.2) + 
  ylim(0,1) + theme(legend.position="top", text = element_text(family = "Palatino", size = 16)) + ylab("Coverage") + xlab("Time") 

data_coverage = data.frame(coverage = c(apply(coverage_fcv, 2, mean), apply(coverage_fcv_online, 2, mean)),
                           time = rep(seq(ncol(coverage_fcv)), 2), Method = rep(c("FCV", "AFCV"), each = ncol(coverage_fcv)))
sd = apply(coverage_fcv, 2, sd)
p2 <- ggplot(data_coverage) + geom_line(aes(x = time, y = coverage, col = Method)) + theme_bw() + geom_hline(yintercept=0.9, linetype="dashed", size=0.8) + 
  geom_ribbon(aes(x = time, ymin = coverage - qnorm(0.95) * sd/sqrt(500), ymax = coverage + qnorm(0.95) * sd/sqrt(500), fill = Method), alpha = 0.2) +
  ylim(0,1) + theme(legend.position="top", text = element_text(family = "Palatino", size = 16)) + ylab("Coverage") + xlab("Time")
ggarrange(p1, p2, ncol=2, nrow=1)


###########################
# Plot time-average coverage
###########################
cummean <- function(x){
  return(cumsum(x)/seq(length(x)))
}
coverage = t(apply(coverage, 1, cummean))
coverage_online = t(apply(coverage_online, 1, cummean))
coverage_fcv = t(apply(coverage_fcv, 1, cummean))
coverage_fcv_online = t(apply(coverage_fcv_online, 1, cummean))

quantile_lo <- function(x){
  return(quantile(x, 0.05))
}
quantile_hi <- function(x){
  return(quantile(x, 0.95))
}


data_coverage = data.frame(coverage = c(apply(coverage, 2, mean), apply(coverage_online, 2, mean)),
                           lo = c(apply(coverage, 2, quantile_lo), apply(coverage_online, 2, quantile_lo)),
                           hi = c(apply(coverage, 2, quantile_hi), apply(coverage_online, 2, quantile_hi)),
                           time = rep(seq(ncol(coverage)), 2), Method = rep(c("QFCV", "AQFCV"), each = ncol(coverage)))
p1 <- ggplot(data_coverage) + 
  geom_hline(yintercept=0.9, linetype="dashed", size=0.8) + 
  geom_ribbon(aes(x = time, y = coverage, ymin = lo, ymax = hi, col = Method, fill = Method), alpha = 0.3) + theme_bw() + 
  theme(legend.position="top", text = element_text(family = "Palatino", size = 16)) + ylab("Coverage") + xlab("Time") 

data_coverage_fcv = data.frame(coverage = c(apply(coverage_fcv, 2, median), apply(coverage_fcv_online, 2, median)),
                              lo = c(apply(coverage_fcv, 2, quantile_lo), apply(coverage_fcv_online, 2, quantile_lo)),
                              hi = c(apply(coverage_fcv, 2, quantile_hi), apply(coverage_fcv_online, 2, quantile_hi)),
                              time = rep(seq(ncol(coverage_fcv)), 2), Method = rep(c("FCV", "AFCV"), each = ncol(coverage_fcv)))
p2 <- ggplot(data_coverage_fcv) + 
  geom_hline(yintercept=0.9, linetype="dashed", size=0.8) + 
  geom_ribbon(aes(x = time, y = coverage, ymin = lo, ymax = hi, col = Method, fill = Method), alpha = 0.3) + theme_bw() + 
  theme(legend.position="top", text = element_text(family = "Palatino", size = 16)) + ylab("Coverage") + xlab("Time")
ggarrange(p1, p2, ncol=2, nrow=1)

###########################
# Run simulation for nonstationary data
###########################
set.seed(100)
n = 3000
n0 = 500 ## number of historical data points before start training 
n_test = 5
p = 20
w = 50
w_train = 40
w_val = 5
lag = 5
lag_fcv = n_test + w_train
gamma = 0.005
gamma_fcv = 15
rep = 500
nsim = floor((n-n0)/lag) - n_test/lag
cc = FALSE
alpha = 0.5
coverage = coverage_online = coverage_fcv = coverage_fcv_online = matrix(NA, nrow = rep, ncol = nsim)
length = length_online = length_fcv = length_fcv_online = matrix(NA, nrow = rep, ncol = nsim)
for(r in 1:rep){
  cat("iter ", r, fill = TRUE)
  data = data_gen_nonstationary(n, p, 0, alpha = alpha)
  err_test = err_qfcv = lo_qfcv = hi_qfcv = lo_qfcv_online = hi_qfcv_online = lo_fcv = hi_fcv = lo_fcv_online = hi_fcv_online = rep(NA, nsim)
  fold_start = floor((n0-w)/lag) 
  err_val_list = err_test_list = c()
  err_fcv_list = c()
  theta = 0
  theta_list = c(0)
  theta_fcv = 0
  err_fcv_list_full = fcv(n, data$x, data$y, w_train+n_test, w_train, lag_fcv)
  for(i in 1:nsim){
    split = n0+(i-1) * lag
    err_test[i] = mean(fit(data$x[(split-w_train+1):split,], data$y[(split-w_train+1):split], data$x[(split+1):(split+n_test),], data$y[(split+1):(split+n_test)]))
    fold_fcv = floor((split-w_train-n_test)/lag_fcv) + 1
    err_fcv_list = err_fcv_list_full[1:fold_fcv]
    err_fcv = mean(err_fcv_list)
    if(cc == FALSE){ 
      se_fcv = sd(err_fcv_list)/sqrt(length(err_fcv_list))
    } 
    else{
      var_fcv = var(err_fcv_list)
      m = floor(length(err_fcv_list) / 50) 
      for(s in 1:m){
        var_fcv = var_fcv + 2 * (1 - s/length(err_fcv_list)) * cov(err_fcv_list[1:(length(err_fcv_list)-s)], err_fcv_list[(s+1):(length(err_fcv_list))])
      }
      var_fcv = (var_fcv < 0) * var(err_fcv_list) + (var_fcv >= 0) * var_fcv
      se_fcv = sqrt(var_fcv/length(err_fcv_list))
    }
    
    lo_fcv[i] = err_fcv - qnorm(0.95) * se_fcv
    hi_fcv[i] = err_fcv + qnorm(0.95) * se_fcv
    lo_fcv_online[i] = err_fcv - (qnorm(0.95) + theta_fcv) * se_fcv 
    hi_fcv_online[i] = err_fcv + (qnorm(0.95) + theta_fcv) * se_fcv 
    theta_fcv = theta_fcv + gamma_fcv * (as.numeric(err_test[i] > hi_fcv_online[i] | err_test[i] < lo_fcv_online[i])-0.1)
    
    result_qfcv = fcv_cal(split, (fold_start + i - 1) * (i > 1), data$x[1:split,], data$y[1:split], w, w_train, w_val, lag, 0, err_val_list, err_test_list)
    result_qfcv_online = fcv_cal(split, (fold_start + i - 1) * (i > 1), data$x[1:split,], data$y[1:split], w, w_train, w_val, lag, theta, err_val_list, err_test_list)
    lo_qfcv_online[i] = result_qfcv_online$err_hat_lo
    hi_qfcv_online[i] = result_qfcv_online$err_hat_hi
    err_val_list = result_qfcv$err_val_list
    err_test_list = result_qfcv$err_test_list
    err_qfcv[i] = result_qfcv$err_hat
    lo_qfcv[i] = result_qfcv$err_hat_lo
    hi_qfcv[i] = result_qfcv$err_hat_hi
    theta = theta + gamma * (as.numeric(err_test[i] > hi_qfcv_online[i] | err_test[i] < lo_qfcv_online[i])-0.1)
    theta_list = c(theta_list, theta)
  }
  coverage[r,] = as.numeric(err_test < hi_qfcv & err_test > lo_qfcv)
  coverage_online[r,] = as.numeric(err_test < hi_qfcv_online & err_test > lo_qfcv_online)
  coverage_fcv[r,] = as.numeric(err_test < hi_fcv & err_test > lo_fcv)
  coverage_fcv_online[r,] = as.numeric(err_test < hi_fcv_online & err_test > lo_fcv_online)
  length[r,] = hi_qfcv - lo_qfcv
  length_online[r,] = hi_qfcv_online - lo_qfcv_online
  length_fcv[r,] = hi_fcv - lo_fcv
  length_fcv_online[r,] = hi_fcv_online - lo_fcv_online
  cat("true test: ", mean(err_test), "coverage: ", mean(coverage[r,]), "coverage_online: ", mean(coverage_online[r,]), "coverage_fcv: ", mean(coverage_fcv[r,]), "coverage_fcv_online: ", mean(coverage_fcv_online[r,]), fill = TRUE)
}

out = list(coverage = coverage, coverage_online = coverage_online, coverage_fcv = coverage_fcv, coverage_fcv_online = coverage_fcv_online,
           length = length, length_online = length_online, length_fcv = length_fcv, length_fcv_online = length_fcv_online)
save(out, file = paste0("AR1_5_nonstationary.RData"))

###########################
# Look at results to plot coverage
###########################
load(file = "AR1_5_nonstationary.RData")
coverage = out$coverage
coverage_online = out$coverage_online
coverage_fcv = out$coverage_fcv
coverage_fcv_online = out$coverage_fcv_online

###########################
# Plot instance-average coverage
###########################
data_coverage = data.frame(coverage = c(apply(coverage, 2, mean), apply(coverage_online, 2, mean)),
                           time = rep(seq(ncol(coverage)), 2), Method = rep(c("QFCV", "AQFCV"), each = ncol(coverage)))
sd = apply(coverage, 2, sd)
p1 <- ggplot(data_coverage) + geom_line(aes(x = time, y = coverage, col = Method)) + theme_bw() + geom_hline(yintercept=0.9, linetype="dashed", size=0.8) + 
  geom_ribbon(aes(x = time, ymin = coverage - qnorm(0.95) * sd/sqrt(500), ymax = coverage + qnorm(0.95) * sd/sqrt(500), fill = Method), alpha = 0.2) + 
  ylim(0,1) + theme(legend.position="top", text = element_text(family = "Palatino", size = 16)) + ylab("Coverage") + xlab("Time") 

data_coverage = data.frame(coverage = c(apply(coverage_fcv, 2, mean), apply(coverage_fcv_online, 2, mean)),
                           time = rep(seq(ncol(coverage_fcv)), 2), Method = rep(c("FCV", "AFCV"), each = ncol(coverage_fcv)))
sd = apply(coverage_fcv, 2, sd)
p2 <- ggplot(data_coverage) + geom_line(aes(x = time, y = coverage, col = Method)) + theme_bw() + geom_hline(yintercept=0.9, linetype="dashed", size=0.8) + 
  geom_ribbon(aes(x = time, ymin = coverage - qnorm(0.95) * sd/sqrt(500), ymax = coverage + qnorm(0.95) * sd/sqrt(500), fill = Method), alpha = 0.2) +
  ylim(0,1) + theme(legend.position="top", text = element_text(family = "Palatino", size = 16)) + ylab("Coverage") + xlab("Time")
ggarrange(p1, p2, ncol=2, nrow=1)


###########################
# Plot time-average coverage
###########################
cummean <- function(x){
  return(cumsum(x)/seq(length(x)))
}
coverage = t(apply(coverage, 1, cummean))
coverage_online = t(apply(coverage_online, 1, cummean))
coverage_fcv = t(apply(coverage_fcv, 1, cummean))
coverage_fcv_online = t(apply(coverage_fcv_online, 1, cummean))

quantile_lo <- function(x){
  return(quantile(x, 0.05))
}
quantile_hi <- function(x){
  return(quantile(x, 0.95))
}


data_coverage = data.frame(coverage = c(apply(coverage, 2, mean), apply(coverage_online, 2, mean)),
                           lo = c(apply(coverage, 2, quantile_lo), apply(coverage_online, 2, quantile_lo)),
                           hi = c(apply(coverage, 2, quantile_hi), apply(coverage_online, 2, quantile_hi)),
                           time = rep(seq(ncol(coverage)), 2), Method = rep(c("QFCV", "AQFCV"), each = ncol(coverage)))
p1 <- ggplot(data_coverage) + 
  geom_hline(yintercept=0.9, linetype="dashed", size=0.8) + 
  geom_ribbon(aes(x = time, y = coverage, ymin = lo, ymax = hi, col = Method, fill = Method), alpha = 0.3) + theme_bw() + 
  theme(legend.position="top", text = element_text(family = "Palatino", size = 16)) + ylab("Coverage") + xlab("Time") 

data_coverage_fcv = data.frame(coverage = c(apply(coverage_fcv, 2, median), apply(coverage_fcv_online, 2, median)),
                               lo = c(apply(coverage_fcv, 2, quantile_lo), apply(coverage_fcv_online, 2, quantile_lo)),
                               hi = c(apply(coverage_fcv, 2, quantile_hi), apply(coverage_fcv_online, 2, quantile_hi)),
                               time = rep(seq(ncol(coverage_fcv)), 2), Method = rep(c("FCV", "AFCV"), each = ncol(coverage_fcv)))
p2 <- ggplot(data_coverage_fcv) + 
  geom_hline(yintercept=0.9, linetype="dashed", size=0.8) + 
  geom_ribbon(aes(x = time, y = coverage, ymin = lo, ymax = hi, col = Method, fill = Method), alpha = 0.3) + theme_bw() + 
  theme(legend.position="top", text = element_text(family = "Palatino", size = 16)) + ylab("Coverage") + xlab("Time")
ggarrange(p1, p2, ncol=2, nrow=1)







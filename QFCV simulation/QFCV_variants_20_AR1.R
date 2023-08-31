##############################################
# Compare different variants of QFCV methods for ARMA(1,0) noise and n_test = 20
##############################################
library(ggplot2)
library(ggpubr)
source("helper.R")

###########################
# Run simulation for n_test = 20 and noise is ARMA(1,0)
###########################
nsim = 500
n = 2000
n_test = 20
p = 20
w = 80
w_train = 40
w_val = 20
lag = 5
alpha_list = c(0.1, 0.3, 0.5, 0.7, 0.9)
for(alpha in alpha_list){
  err_test = err_fcv_cal = lo_fcv_cal = hi_fcv_cal = err_fcv_cal_mid = lo_fcv_cal_mid = hi_fcv_cal_mid = err_fcv_cal_full = 
    lo_fcv_cal_full = hi_fcv_cal_full = err_fcv_cal_0 = lo_fcv_cal_0 = hi_fcv_cal_0 = err_nfcv = se_nfcv = err_fcv = se_fcv = rep(NA, nsim)
  for(i in 1:nsim){
    data = data_gen(n, p, n_test, alpha = alpha)
    err_test[i] = mean(fit(data$x[(n-w_train+1):n,], data$y[(n-w_train+1):n], data$x_test, data$y_test))
    result_fcv_cal = fcv_cal(n, data$x, data$y, w, w_train, w_val, lag)
    err_fcv_cal[i] = result_fcv_cal$err_hat
    lo_fcv_cal[i] = result_fcv_cal$err_hat_lo
    hi_fcv_cal[i] = result_fcv_cal$err_hat_hi
    err_fcv_cal_mid[i] = result_fcv_cal$err_hat_mid
    lo_fcv_cal_mid[i] = result_fcv_cal$err_hat_mid_lo
    hi_fcv_cal_mid[i] = result_fcv_cal$err_hat_mid_hi
    err_fcv_cal_full[i] = result_fcv_cal$err_hat_full
    lo_fcv_cal_full[i] = result_fcv_cal$err_hat_full_lo
    hi_fcv_cal_full[i] = result_fcv_cal$err_hat_full_hi
    err_fcv_cal_0[i] = result_fcv_cal$err_hat_0
    lo_fcv_cal_0[i] = result_fcv_cal$err_hat_0_lo
    hi_fcv_cal_0[i] = result_fcv_cal$err_hat_0_hi
    
    err_fcv_list = fcv(n, data$x, data$y, w_train+n_test, w_train, w_train+n_test)
    err_fcv[i] = mean(err_fcv_list)
    se_fcv[i] = sd(err_fcv_list)/sqrt(length(err_fcv_list))
    cat("true test: ", err_test[i], "quantile: ", c(lo_fcv_cal[i], hi_fcv_cal[i]), fill = TRUE)
  }
  
  out = list(err_test = err_test, err_fcv = err_fcv, se_fcv = se_fcv,
             err_fcv_cal = err_fcv_cal, lo_fcv_cal = lo_fcv_cal, hi_fcv_cal = hi_fcv_cal,
             err_fcv_cal_0 = err_fcv_cal_0, lo_fcv_cal_0 = lo_fcv_cal_0, hi_fcv_cal_0 = hi_fcv_cal_0,
             err_fcv_cal_mid = err_fcv_cal_mid, lo_fcv_cal_mid = lo_fcv_cal_mid, hi_fcv_cal_mid = hi_fcv_cal_mid,
             err_fcv_cal_full = err_fcv_cal_full, lo_fcv_cal_full = lo_fcv_cal_full, hi_fcv_cal_full = hi_fcv_cal_full)
  save(out, file = paste0("sparse_linear_20_AR1_alpha-", alpha,".RData"))
  print(paste0("Results for alpha-", alpha, " saved to disk."))
}

###########################
# Plot for n_test = 20 and noise is ARMA(1,0)
###########################
alpha_list = c(0.1,0.3,0.5,0.7,0.9)
coverage_cal = c()
length_cal = c()
coverage_cal_0 = c()
length_cal_0 = c()
coverage_cal_mid = c()
length_cal_mid = c()
coverage_cal_full = c()
length_cal_full = c()
lo_length_cal = hi_length_cal = lo_length_cal_0 = hi_length_cal_0 = lo_length_cal_mid = hi_length_cal_mid = lo_length_cal_full = hi_length_cal_full =c()
se_length_cal = se_length_cal_0 = se_length_cal_mid = se_length_cal_full = c()
for(alpha in alpha_list){
  load(paste0("sparse_linear_20_AR1_alpha-", alpha,".RData"))
  err_test = out$err_test
  err_fcv_cal = out$err_fcv_cal
  lo_fcv_cal = out$lo_fcv_cal
  hi_fcv_cal = out$hi_fcv_cal
  err_fcv_cal_0 = out$err_fcv_cal_0
  lo_fcv_cal_0 = out$lo_fcv_cal_0
  hi_fcv_cal_0 = out$hi_fcv_cal_0
  err_fcv_cal_mid = out$err_fcv_cal_mid
  lo_fcv_cal_mid = out$lo_fcv_cal_mid
  hi_fcv_cal_mid = out$hi_fcv_cal_mid
  err_fcv_cal_full = out$err_fcv_cal_full
  lo_fcv_cal_full = out$lo_fcv_cal_full
  hi_fcv_cal_full = out$hi_fcv_cal_full
  length_base = quantile(err_test, 0.95) - quantile(err_test, 0.05)
  coverage_cal <- c(coverage_cal, mean((err_test > lo_fcv_cal) & (err_test < hi_fcv_cal)))
  length_cal <- c(length_cal, mean(hi_fcv_cal - lo_fcv_cal)/length_base)
  coverage_cal_0 <- c(coverage_cal_0, mean((err_test > lo_fcv_cal_0) & (err_test < hi_fcv_cal_0)))
  length_cal_0 <- c(length_cal_0, mean(hi_fcv_cal_0 - lo_fcv_cal_0)/length_base)
  coverage_cal_mid <- c(coverage_cal_mid, mean((err_test > lo_fcv_cal_mid) & (err_test < hi_fcv_cal_mid)))
  length_cal_mid <- c(length_cal_mid, mean(hi_fcv_cal_mid - lo_fcv_cal_mid)/length_base)
  coverage_cal_full <- c(coverage_cal_full, mean((err_test > lo_fcv_cal_full) & (err_test < hi_fcv_cal_full)))
  length_cal_full <- c(length_cal_full, mean(hi_fcv_cal_full - lo_fcv_cal_full)/length_base)
  se_length_cal <- as.vector(c(se_length_cal, sd(hi_fcv_cal - lo_fcv_cal)/length_base/sqrt(nsim)))
  se_length_cal_0 <- as.vector(c(se_length_cal_0, sd(hi_fcv_cal_0 - lo_fcv_cal_0)/length_base/sqrt(nsim)))
  se_length_cal_mid <- as.vector(c(se_length_cal_mid, sd(hi_fcv_cal_mid - lo_fcv_cal_mid)/length_base/sqrt(nsim)))
  se_length_cal_full <- as.vector(c(se_length_cal_full, sd(hi_fcv_cal_full - lo_fcv_cal_full)/length_base/sqrt(nsim)))
}

data <- data.frame(coverage = c(coverage_cal, coverage_cal_0, coverage_cal_mid, coverage_cal_full), alpha = rep(alpha_list, 4),
                   length = c(length_cal, length_cal_0, length_cal_mid, length_cal_full),
                   Method = rep(c("QFCV(1)", "QFCV(0)", "QFCV(2)", "QFCV(20)"), each = length(coverage_cal)))
p1 <- ggplot(data) + geom_line(aes(x = alpha, y = coverage, col = Method), size = 0.8) + 
  geom_ribbon(aes(ymin = coverage - qnorm(0.95) * sqrt(coverage * (1-coverage)/nsim), 
                  ymax = coverage + qnorm(0.95) * sqrt(coverage * (1-coverage)/nsim),
                  x = alpha, fill = Method), alpha = 0.1) +
  theme_bw() + ylim(c(0.7,1)) + xlab("Autocorrelation") + ylab("Coverage") + theme(text = element_text(family = "Palatino", size = 16))
p2 <- ggplot(data) + geom_line(aes(x = alpha, y = length, col = Method), size = 0.8) +
  geom_ribbon(aes(ymin = length - qnorm(0.95) * c(se_length_cal, se_length_cal_0, se_length_cal_mid, se_length_cal_full),
                  ymax = length + qnorm(0.95) * c(se_length_cal, se_length_cal_0, se_length_cal_mid, se_length_cal_full),
                  x = alpha, fill = Method), alpha = 0.1) + ylim(c(0.6, 1.1)) + 
  theme_bw() + xlab("Autocorrelation") + ylab("Length Ratio") + theme(text = element_text(family = "Palatino", size = 16))
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE)

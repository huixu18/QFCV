##############################################
# Compare different variants of FCV methods for ARMA(1,20) noise and n_test = 20
##############################################
library(ggplot2)
library(ggpubr)
source("helper.R")

###########################
# Run simulation for n_test = 20 and noise is ARMA(1,20)
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
  ## err_test: true test error 
  ## err_fcv_cal: point estimator for QFCV 
  ## [lo_fcv_cal, hi_fcv_cal]: predictive interval for QFCV 
  ## err_fcv: point estimator for FCV 
  ## se_fcv, se_fcv_cc, se_fcv_pi, se_fcv_pi_cc: se for naive FCV, FCV with covariance correction, FCV as predictive interval, FCV with covariance correction as predictive interval
  err_test = err_fcv_cal = lo_fcv_cal = hi_fcv_cal = err_fcv = se_fcv = se_fcv_cc = se_fcv_pi = se_fcv_pi_cc = rep(NA, nsim)
  for(i in 1:nsim){
    data = data_gen(n, p, n_test, alpha = alpha, MA = TRUE)
    err_test[i] = mean(fit(data$x[(n-w_train+1):n,], data$y[(n-w_train+1):n], data$x_test, data$y_test))
    
    ## QFCV predictive interval
    result_fcv_cal = fcv_cal(n, data$x, data$y, w, w_train, w_val, lag)
    err_fcv_cal[i] = result_fcv_cal$err_hat
    lo_fcv_cal[i] = result_fcv_cal$err_hat_lo
    hi_fcv_cal[i] = result_fcv_cal$err_hat_hi
    
    ## naive FCV confidence interval and predictive interval
    err_fcv_list = fcv(n, data$x, data$y, w_train+n_test, w_train, n_test)
    err_fcv[i] = mean(err_fcv_list)
    se_fcv[i] = sd(err_fcv_list)/sqrt(length(err_fcv_list))
    se_fcv_pi[i] = sd(err_fcv_list)
    
    ## FCV confidence and predictive interval with covariance correction
    var_fcv = var(err_fcv_list)
    m = floor(length(err_fcv_list) / 50) 
    for(s in 1:m){
      var_fcv = var_fcv + 2 * (1 - s/length(err_fcv_list)) * cov(err_fcv_list[1:(length(err_fcv_list)-s)], err_fcv_list[(s+1):(length(err_fcv_list))])
    }
    if(var_fcv < 0){
      var_fcv = var(err_fcv_list)
    }
    se_fcv_cc[i] = sqrt(var_fcv/length(err_fcv_list))
    se_fcv_pi_cc[i] = sqrt(var_fcv)
    
    cat("true test: ", err_test[i], "quantile: ", c(lo_fcv_cal[i], hi_fcv_cal[i]), fill = TRUE)
  }
  
  out = list(err_test = err_test, err_fcv = err_fcv, se_fcv = se_fcv,
             err_fcv_cal = err_fcv_cal, lo_fcv_cal = lo_fcv_cal, hi_fcv_cal = hi_fcv_cal,
             err_fcv = err_fcv, se_fcv = se_fcv, se_fcv_cc = se_fcv_cc, se_fcv_pi = se_fcv_pi, se_fcv_pi_cc = se_fcv_pi_cc)
  save(out, file = paste0("sparse_linear_5_AR1_alpha-", alpha,".RData"))
  print(paste0("Results for alpha-", alpha, " saved to disk."))
}

###########################
# Plot for n_test = 20 and noise is ARMA(1,20)
###########################
alpha_list = c(0.1,0.3,0.5,0.7,0.9)
coverage_cal = c()
length_cal = c()
coverage_fcv = c()
length_fcv = c()
coverage_fcv_pi = c()
length_fcv_pi = c()
coverage_fcv_cc = c()
length_fcv_cc = c()
coverage_fcv_pi_cc = c()
length_fcv_pi_cc = c()

lo_length_cal = hi_length_cal = se_length_cal = se_length_fcv = se_length_fcv_cc = se_length_fcv_pi = se_length_fcv_pi_cc = c()
for(alpha in alpha_list){
  load(paste0("sparse_linear_20_AR1_MA20_alpha-", alpha,".RData"))
  err_test = out$err_test
  err_fcv_cal = out$err_fcv_cal
  lo_fcv_cal = out$lo_fcv_cal
  hi_fcv_cal = out$hi_fcv_cal
  err_fcv = out$err_fcv
  se_fcv = out$se_fcv
  se_fcv_cc = out$se_fcv_cc 
  se_fcv_pi = out$se_fcv_pi
  se_fcv_pi_cc = out$se_fcv_pi_cc
  length_base = quantile(err_test, 0.95) - quantile(err_test, 0.05)
  coverage_cal <- c(coverage_cal, mean((err_test > lo_fcv_cal) & (err_test < hi_fcv_cal)))
  length_cal <- c(length_cal, mean(hi_fcv_cal - lo_fcv_cal)/length_base)
  coverage_fcv = c(coverage_fcv, mean((err_test > err_fcv - qnorm(0.95) * se_fcv) & (err_test < err_fcv + qnorm(0.95) * se_fcv)))
  length_fcv <- c(length_fcv, 2 * qnorm(0.95) * mean(se_fcv)/length_base)
  coverage_fcv_cc = c(coverage_fcv_cc, mean((err_test > err_fcv - qnorm(0.95) * se_fcv_cc) & (err_test < err_fcv + qnorm(0.95) * se_fcv_cc)))
  length_fcv_cc = c(length_fcv_cc, 2 * qnorm(0.95) * mean(se_fcv_cc)/length_base)
  coverage_fcv_pi = c(coverage_fcv_pi, mean((err_test > err_fcv - qnorm(0.95) * se_fcv_pi) & (err_test < err_fcv + qnorm(0.95) * se_fcv_pi)))
  length_fcv_pi = c(length_fcv_pi, 2 * qnorm(0.95) * mean(se_fcv_pi)/length_base)
  coverage_fcv_pi_cc = c(coverage_fcv_pi_cc, mean((err_test > err_fcv - qnorm(0.95) * se_fcv_pi_cc) & (err_test < err_fcv + qnorm(0.95) * se_fcv_pi_cc)))
  length_fcv_pi_cc = c(length_fcv_pi_cc, 2 * qnorm(0.95) * mean(se_fcv_pi_cc)/length_base)
  
  se_length_cal <- as.vector(c(se_length_cal, sd(hi_fcv_cal - lo_fcv_cal)/length_base/sqrt(500)))
  se_length_fcv <- as.vector(c(se_length_fcv, sd(2 * qnorm(0.95) * se_fcv)/length_base/sqrt(500)))
  se_length_fcv_cc <- as.vector(c(se_length_fcv_cc, sd(2 * qnorm(0.95) * se_fcv_cc)/length_base/sqrt(500)))
  se_length_fcv_pi <- as.vector(c(se_length_fcv_pi, sd(2 * qnorm(0.95) * se_fcv_pi)/length_base/sqrt(500)))
  se_length_fcv_pi_cc <- as.vector(c(se_length_fcv_pi_cc, sd(2 * qnorm(0.95) * se_fcv_pi_cc)/length_base/sqrt(500)))
}

data <- data.frame(coverage = c(coverage_cal, coverage_fcv, coverage_fcv_cc, coverage_fcv_pi), alpha = rep(alpha_list, 4),
                   length = c(length_cal, length_fcv, length_fcv_cc, length_fcv_pi),
                   Method = rep(c("QFCV", "FCV", "FCV(c)", "FCV(p)"), each = length(coverage_cal)))
p1 <- ggplot(data) + geom_line(aes(x = alpha, y = coverage, col = Method), size = 0.8) + geom_hline(yintercept=0.9, linetype="dashed", size = 1) + 
  geom_ribbon(aes(ymin = coverage - qnorm(0.95) * sqrt(coverage * (1-coverage)/500), 
                  ymax = coverage + qnorm(0.95) * sqrt(coverage * (1-coverage)/500),
                  x = alpha, fill = Method), alpha = 0.2) +
  theme_bw() + xlab("Autocorrelation") + ylab("Coverage") + theme(text = element_text(family = "Palatino", size = 16))
p2 <- ggplot(data) + geom_line(aes(x = alpha, y = length, col = Method), size = 0.8) +
  geom_ribbon(aes(ymin = length - qnorm(0.95) * c(se_length_cal, se_length_fcv, se_length_fcv_cc, se_length_fcv_pi),
                  ymax = length + qnorm(0.95) * c(se_length_cal, se_length_fcv, se_length_fcv_cc, se_length_fcv_pi),
                  x = alpha, fill = Method), alpha = 0.2) +
  theme_bw() + xlab("Autocorrelation") + ylab("Length Ratio") + theme(text = element_text(family = "Palatino", size = 16))
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE)











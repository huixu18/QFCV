##############################################
# Plot comparison of stochastic and expected test errors, along with naive FCV and QFCV intervals for 50 simulation instances.
##############################################
library(ggplot2)
source("helper.R")


##############################################
# Run the simulation 
##############################################
nsim = 500
n = 1000
n_test = 20
p = 20
w = 80
w_train = 40
w_val = 20
lag = 5
alpha = 0.5
err_test = err_fcv_cal = lo_fcv_cal = hi_fcv_cal = err_fcv = se_fcv = se_fcv_cc = se_fcv_pi = rep(NA, nsim)
for(i in 1:nsim){
  data = data_gen(n, p, n_test, alpha = alpha)
  err_test[i] = mean(fit(data$x[(n-w_train+1):n,], data$y[(n-w_train+1):n], data$x_test, data$y_test))
  result_fcv_cal = fcv_cal(n, data$x, data$y, w, w_train, w_val, lag)
  err_fcv_cal[i] = result_fcv_cal$err_hat
  lo_fcv_cal[i] = result_fcv_cal$err_hat_lo
  hi_fcv_cal[i] = result_fcv_cal$err_hat_hi
  
  err_fcv_list = fcv(n, data$x, data$y, w_train+n_test, w_train, n_test)
  err_fcv[i] = mean(err_fcv_list)
  se_fcv[i] = sd(err_fcv_list)/sqrt(length(err_fcv_list))
  se_fcv_pi[i] = sd(err_fcv_list)
  var_fcv = var(err_fcv_list)
  m = floor(length(err_fcv_list) / 20) 
  for(s in 1:m){
    var_fcv = var_fcv + 2 * (1 - s/length(err_fcv_list)) * cov(err_fcv_list[1:(length(err_fcv_list)-s)], err_fcv_list[(s+1):(length(err_fcv_list))])
  }
  if(var_fcv < 0){
    var_fcv = var(err_fcv_list)
  }
  se_fcv_cc[i] = sqrt(var_fcv/length(err_fcv_list))
  cat("true test: ", err_test[i], "quantile: ", c(lo_fcv_cal[i], hi_fcv_cal[i]), "se_fcv: ", se_fcv[i], "m: ", m, fill = TRUE)
}

out = list(err_test = err_test, err_fcv = err_fcv, se_fcv = se_fcv, se_fcv_cc = se_fcv_cc, se_fcv_pi = se_fcv_pi,
           err_fcv_cal = err_fcv_cal, lo_fcv_cal = lo_fcv_cal, hi_fcv_cal = hi_fcv_cal)
save(out, file = paste0("sparse_linear_20_alpha-", alpha,".RData"))
print(paste0("Results for alpha-", alpha, " saved to disk."))

##############################################
# Look at result and generate plot
##############################################
load(paste0("sparse_linear_20_alpha-", alpha,".RData"))
v = 50
err_test = out$err_test[1:v]
err_fcv_cal = out$err_fcv_cal[1:v]
lo_fcv_cal = out$lo_fcv_cal[1:v]
hi_fcv_cal = out$hi_fcv_cal[1:v]
err_fcv = out$err_fcv[1:v]
se_fcv = out$se_fcv[1:v] 
se_fcv_cc = out$se_fcv_cc[1:v]
se_fcv_pi = out$se_fcv_pi[1:v]

data_PI = data.frame(True = rep(err_test, 2), lo = c(err_fcv - qnorm(0.95) * se_fcv, lo_fcv_cal), hi = c(err_fcv + qnorm(0.95)* se_fcv, hi_fcv_cal),
                     Method = rep(c("FCV","QFCV"), each = length(err_test)), index = rep(seq(length(err_test)), 2),
                     estimate = rep(err_fcv, 2))
ggplot(data_PI) + geom_point(aes(x = index, y = True)) + geom_segment(aes(x = index, xend = index, y = lo, yend = hi, col = Method), size = 2, alpha = 0.4) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + 
  theme_bw() + ylab("True Test Error") + xlab("Index") + geom_hline(yintercept=mean(err_test), linetype="dashed", size=1) + 
  geom_hline(yintercept = quantile(err_test, 0.05), linetype = "dashed", size = 1, col = "orange") + 
  geom_hline(yintercept = quantile(err_test, 0.95), linetype = "dashed", size = 1, col = "orange") + theme(text = element_text(family = "Palatino", size = 16)) 

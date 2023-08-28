library(ggplot2)
source("helper.R")

set.seed(100)
nsim = 1000
n = 40
n_test = 5
p = 20
w = 50
w_train = 40
lag = 5
err_test = rep(NA, nsim)
for(i in 1:nsim){
  data = data_gen(n, p, n_test, alpha = 0.5)
  err_test[i] = mean(fit(data$x[(n-w_train+1):n,], data$y[(n-w_train+1):n], data$x_test, data$y_test))
}
data = data_gen(2000, p, n_test, alpha = 0.5)
err_fcv_list = fcv(2000, data$x, data$y, w_train+n_test, w_train, lag)
nwindow = length(err_fcv_list)

result = data.frame(Err = c(err_test, err_fcv_list[1:30], err_fcv_list[1:300]), Type = c(rep("True", nsim), rep("K = 30", 30), rep("K = 300", 300)))
ggplot(result, aes(x = Err, fill = Type)) + geom_histogram(aes(y = ..density..), binwidth = 0.5, alpha = 0.6) + 
  geom_density(aes(col = Type), size = 1, alpha = 0.1) + ylab("Density") + xlab("") + 
  theme_bw() + theme(axis.text.y=element_blank()) + facet_wrap(~Type,nrow=3) + theme(text = element_text(family = "Palatino", size = 16)) 

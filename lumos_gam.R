library(mboost)
library(gbm)
library(ggplot2)

#Run script to load and clean data
source("lumos_cleaning.R")

##### Fit initial BAMMO model #####
gammo_form = Throughput ~ bbs(movingSpeed) + bbs(compassDirection) + 
                            bspatial(latitude, longitude) + 
                            bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) + 
                            bols(nrStatus, contrasts.arg = "contr.sum") +
                            bbs(lte_rssi) +  bbs(lte_rsrp) + bbs(lte_rsrq) +
                            bbs(nr_ssRsrp) + bbs(nr_ssSinr) 

gammo = gamboost(gammo_form,
                 data = lumos, 
                 na.action = na.pass,
                 control = boost_control(mstop = 100, nu = 0.1)) 

# Tune number of boosting iterations with 10-fold CV
set.seed(2023)
cv_gammo = cvrisk(gammo, folds = cv(weights = model.weights(gammo), type = "kfold", B = 10))

# Find smallest mstop with average risk within 1 se of min risk
cm = apply(cv_gammo, 2, mean)
cse = apply(cv_gammo, 2, sd)/sqrt(nrow(cv_gammo))
mstops = attr(cv_gammo, "mstop")
mstop_min = mstops[which.min(cm)]
ind_1se = which(cm < min(cm) + cse[which.min(cm)])[1]
mstop_1se = mstops[ind_1se]

c(mstop_min, min(cm))
c(mstop_1se, cm[ind_1se])

#Plot CV curve and save to pdf
plot(cv_gammo, main = "", xlab = "Number of boosting iterations (M)", ylab = "Mean Squared Error (MSE)")
lines(c(ind_1se, ind_1se), c(0, 2*min(cm)), lty = 1, lwd = 1)
lines(c(mstop_min+1, mstop_min+1), c(0, 2*min(cm)), lty = 2, lwd = 1)

pdf("gam_M_cv.pdf", width = 7, height = 5)
plot(cv_gammo, main = "", xlab = "Number of boosting iterations (M)", ylab = "Mean Squared Error (MSE)")
lines(c(ind_1se, ind_1se), c(0, 2*min(cm)), lty = 1, lwd = 1)
lines(c(mstop_min+1, mstop_min+1), c(0, 2*min(cm)), lty = 2, lwd = 1)
dev.off()

#variable importance (MSE reduction) plot for M_min model
gammo_varimp = varimp(gammo)
pos_varimp = sort(gammo_varimp[gammo_varimp > 0],decreasing = F)
varimp_per = pos_varimp/sum(abs(pos_varimp))
varimp_names = c("5G NR status", "(latitude, longitude)", "moving speed", "throughput (lag 2)", "5G NR RSRP", "throughput (lag 3)", "5G NR SSNR","throughput (lag 1)")

ggplot(mapping = aes(x = varimp_names, y = varimp_per*100)) + 
    geom_col() + 
    scale_x_discrete(limits = varimp_names) + 
    scale_y_continuous(limits = c(0,100), expand = c(0,0)) + 
    coord_flip() + 
    theme_bw() +
    theme(plot.margin = margin(10, 20, 10, 10, "points")) +
    labs(y = "MSE Reduction (%)", x = "Variable")

ggsave("gam_varimp.pdf", width = 7, height = 3)

# count variable selections in M_min model
names(gammo$coef())
gammo$coef()
gammo$xselect()

x_selections = names(gammo$baselearner)[gammo$xselect()]
var_list = unique(x_selections)
var_count = rep(0, length(var_list))
for(i in 1:length(var_list)){
    var_count[i] = sum(x_selections == var_list[i])
}
cbind(var_list, var_count)

# Update M for BAMMO (using 1 se rule)
gammo[mstop_1se]

# Count variable selections in M_1se model
names(gammo$coef())
gammo$coef()
gammo$xselect()

x_selections = names(gammo$baselearner)[gammo$xselect()]
var_list = unique(x_selections)
var_count = rep(0, length(var_list))
for(i in 1:length(var_list)){
    var_count[i] = sum(x_selections == var_list[i])
}
cbind(var_list, var_count)

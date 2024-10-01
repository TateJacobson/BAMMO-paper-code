library(mboost)

#Run script to load and clean data
source("lumos_cleaning.R")

gammo_form = Throughput ~ bbs(movingSpeed) + 
    bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) + 
    bbs(nr_ssRsrp) + bbs(nr_ssSinr) 

gammo = gamboost(gammo_form,
                 data = lumos, 
                 na.action = na.pass,
                 control = boost_control(mstop = 100, nu = 0.1))

# Tune number of boosting iterations with 10-fold CV
set.seed(2023)
cv_gammo = cvrisk(gammo, folds = cv(weights = model.weights(gammo), type = "kfold", B = 10))
cm = apply(cv_gammo, 2, mean)
mstops = attr(cv_gammo, "mstop")
mstop_min = mstops[which.min(cm)]
gammo[mstop_min]

# plot estimated partial effects #####
names(gammo$baselearner) #use indices for which below

# get estimated offsets for NA values
rsrp_na_pred <- predict(gammo, newdata = data.frame(nr_ssRsrp = c( mean( lumos$nr_ssRsrp, na.rm = T)  , NA) ), which = 5)[2] #5G NR RSRP
ssnr_na_pred <- predict(gammo, newdata = data.frame(nr_ssSinr = c( mean( lumos$nr_ssSinr, na.rm = T ), NA) ), which = 6)[2] #5G NR SSNR

# plot effects side-by-side and with same y limits
pdf("partial_effects.pdf", width = 8, height = 4)

par(mfrow = c(1, 3), mar = c(4, 5, 0.1, 0.1))
plot(gammo, which = 5, rug = F, newdata = data.frame(nr_ssRsrp = seq( min(lumos$nr_ssRsrp, na.rm = T), max(lumos$nr_ssRsrp, na.rm = T), length = 20 ) ), 
     type = "b", xlab = "5G NR RSRP", ylim = c(-55,105))
abline(h = rsrp_na_pred, lty = 2)
plot(gammo, which = 6, rug = F, newdata = data.frame(nr_ssSinr = seq( min(lumos$nr_ssSinr, na.rm = T), max(lumos$nr_ssSinr, na.rm = T), length = 20 ) ),
     type = "b", xlab = "5G NR RSSNR", ylim = c(-55,105))
abline(h = ssnr_na_pred, lty = 2)
plot(gammo, which = 1, rug = F, newdata = data.frame(movingSpeed = seq( min(lumos$movingSpeed, na.rm = T), max(lumos$movingSpeed, na.rm = T), length = 20 ) ),
     type = "b", xlab = "Moving Speed", ylim = c(-55,105))
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))

dev.off()


#plot partial effects for lags
pdf("lag_partial_effects.pdf", width = 8, height = 4)

par(mfrow = c(1, 3), mar = c(4, 5, 0.1, 0.1))
plot(gammo, which = 2, rug = F, 
     newdata = data.frame(Throughput_lag1 = seq( min(lumos$Throughput, na.rm = T), max(lumos$Throughput, na.rm = T), length = 20 ) ), 
     type = "b", xlab = "Throughput (lag 1)")
plot(gammo, which = 3, rug = F, 
     newdata = data.frame(Throughput_lag2 = seq( min(lumos$Throughput, na.rm = T), max(lumos$Throughput, na.rm = T), length = 20 ) ), 
     type = "b", xlab = "Throughput (lag 2)")
plot(gammo, which = 4, rug = F, 
     newdata = data.frame(Throughput_lag3 = seq( min(lumos$Throughput, na.rm = T), max(lumos$Throughput, na.rm = T), length = 20 ) ), 
     type = "b", xlab = "Throughput (lag 3)")
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))

dev.off()

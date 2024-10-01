library(mice)
library(missForest)
library(missRanger)
library(VIM)
library(microbenchmark)
library(mboost)
library(gbm)

##### Run script to load and clean data #####
source("lumos_cleaning.R")

rm_vars <- c("run_num", "seq_num", "lte_rssnr", "nr_ssRsrq", "Throughput" )
X <- lumos[, !(colnames(lumos) %in% rm_vars)]

Y <- lumos$Throughput

#get a subset of the data for the timing comparison
n <-nrow(X)
samp <- sample(1:n, floor(n/10))
Y <- Y[samp]
X <- X[samp,]

n <- nrow(X)

#split into training and test data
train_ind <- sample(1:n, floor(n*3/4))

X_train <- X[train_ind,]
Y_train <- Y[train_ind]
X_test <- X[-train_ind,]
Y_test <- Y[-train_ind]

##### Compare timing of methods of handling missing data #####
mstop_init <- 100
inner_folds <- 10

#mean imputation
imp_time <- microbenchmark({
    imp_X_train <- X_train
    imp_X_test <- X_test
    
    for(i in 1:ncol(imp_X_train)){
        if(is.numeric(imp_X_train[,i])){
            var_mean <- mean(imp_X_train[,i], na.rm = T)
        } else if (is.factor(imp_X_train[,i])){
            var_mean <- names(which.max(table(imp_X_train[,i])))
        }
        imp_X_train[is.na(imp_X_train[,i]), i] <- var_mean
        imp_X_test[is.na(imp_X_test[,i]), i] <- var_mean
    }
    
    imp_gam_train <- gamboost(Y_train ~ bbs(movingSpeed) + bbs(compassDirection) + #M
                                 bspatial(latitude, longitude) + #L
                                 bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) + 
                                 bbs(lte_rssi) +  bbs(lte_rsrp) + bbs(lte_rsrq) +
                                 bbs(nr_ssRsrp) + bbs(nr_ssSinr), 
                             data = imp_X_train, 
                             na.action = na.omit,
                             control = boost_control(mstop = mstop_init))
    
    cv_imp_gam <- cvrisk(imp_gam_train, folds = cv(weights = model.weights(imp_gam_train), type = "kfold", B = inner_folds))
    imp_gam_train[mstop(cv_imp_gam)]
    imp_gam_preds <- imp_gam_train$predict(imp_X_test)
}, times = 10, unit = "second") 

#MICE
mice_imps <- 5

mice_time <- microbenchmark({
    
    mice_train <- mice(X_train, printFlag = F, m = mice_imps)
    mice_test <- mice(X_test, printFlag = F, m = mice_imps)
    
    mice_pred_mat <- matrix(NA, nrow = nrow(X_test), ncol = mice_imps)
    
    for(m in 1:mice_imps){
        mice_train_m <- data.frame(Throughput = Y_train, complete(mice_train, m))
        
        mice_gam_train <- gamboost(Throughput ~ bbs(movingSpeed) + bbs(compassDirection) + #M
                                       bspatial(latitude, longitude) + #L
                                       bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) + 
                                       bbs(lte_rssi) +  bbs(lte_rsrp) + bbs(lte_rsrq) +
                                       bbs(nr_ssRsrp) + bbs(nr_ssSinr),
                                   data = mice_train_m,
                                   na.action = na.omit,
                                   control = boost_control(mstop = mstop_init))
        
        cv_mice_gam <- cvrisk(mice_gam_train, folds = cv(weights = model.weights(mice_gam_train), type = "kfold", B = inner_folds))
        mice_gam_train[mstop(cv_mice_gam)]
        mice_pred_mat[, m] <- mice_gam_train$predict(complete(mice_test, m)) 
    }
}, times = 10, unit = "second") 

#missForest 
ranger_time <- microbenchmark({
    lumos_missRanger_train <- missRanger(X_train, num.trees = 50, data_only = T, verbose = F)
    
    missRanger_gam_train <- gamboost(Y_train ~ bbs(movingSpeed) + bbs(compassDirection) + #M
                                         bspatial(latitude, longitude) + #L
                                         bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) + 
                                         bbs(lte_rssi) +  bbs(lte_rsrp) + bbs(lte_rsrq) +
                                         bbs(nr_ssRsrp) + bbs(nr_ssSinr),
                                     data = lumos_missRanger_train,
                                     na.action = na.omit,
                                     control = boost_control(mstop = mstop_init))
    
    cv_missRanger_gam <- cvrisk(missRanger_gam_train, folds = cv(weights = model.weights(missRanger_gam_train), type = "kfold", B = inner_folds))
    missRanger_gam_train[mstop(cv_missRanger_gam)]
    
    lumos_missRanger_test <- missRanger(X_test, num.trees = 50, data_only = T, verbose = F)
    missRanger_gam_preds <- missRanger_gam_train$predict(lumos_missRanger_test)
}, times = 10, unit = "second") 

# BAMMO
gammo_time <- microbenchmark({
    gam_train <- gamboost(Y_train ~ bbs(movingSpeed) + bbs(compassDirection) +
                         bspatial(latitude, longitude) + 
                         bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) + 
                         bols(nrStatus, contrasts.arg = "contr.sum") +
                         bbs(lte_rssi) +  bbs(lte_rsrp) + bbs(lte_rsrq) +
                         bbs(nr_ssRsrp) + bbs(nr_ssSinr), 
                     data = X_train, 
                     na.action = na.pass,
                     control = boost_control(mstop = mstop_init))

    cv_gam_train <- cvrisk(gam_train, folds = cv(weights = model.weights(gam_train), type = "kfold", B = inner_folds))
    gam_train[mstop(cv_gam_train)]
    gam_preds <- gam_train$predict(X_test)
}, times = 10, unit = "second")

##### Make a table of times #####
time_table <- data.frame(BAMMO = gammo_time$time, IMP = imp_time$time, MICE = mice_time$time, missForest = ranger_time$time)
time_table <- time_table/1e9
time_table

time_means <- apply(time_table, 2, mean)
time_ses <- apply(time_table, 2, function(c) sd(c)/sqrt( length(c) )  )

library(mboost)
library(gbm)
library(mice)
library(missRanger)
library(parallel)

ppn <- as.numeric(Sys.getenv("SLURM_NTASKS"))
ncores <- ifelse(is.na(ppn), 1, ppn)

#Run script to load and clean data
source('lumos_cleaning.R')

#### Prediction comparison with CV ####
lumos_cv = function(df, reps = 100, test_prop = 1/10, inner_folds = 10, mstop_init = 100, 
                    train_na_vars = NULL, test_na_vars = NULL, train_na_prop = 0, test_na_prop = 0,
                    mice_imps = 5){
    
    # Variables to remove from imputation
    rm_vars <- c("run_num", "seq_num", "lte_rssnr", "nr_ssRsrq", "Throughput" )
    
    sim_loop = function(r){

        tryagain = FALSE
        try_iter = 1

        repeat {
            tryagain = tryCatch({
                # pull training and test data
                n <- nrow(df)
                test_inds <- sample(n, n*test_prop)
                
                test_df <- df[test_inds,]
                train_df <- df[-test_inds,]
                
                # Censor observations in training and test data
                test_inds = which( complete.cases(test_df[, colnames(test_df) %in% test_na_vars]) )
                train_inds = which( complete.cases(train_df[, colnames(train_df) %in% train_na_vars]) )
                
                test_na_inds = sample(test_inds, size = nrow(test_df)*test_na_prop)
                train_na_inds = sample(train_inds, size = nrow(train_df)*train_na_prop)
                
                test_df[test_na_inds, colnames(test_df) %in% test_na_vars] = NA
                train_df[train_na_inds, colnames(train_df) %in% train_na_vars] = NA
                
                # Fit models to training data and tune mstop
                ## BAMMO
                gam_train = gamboost(Throughput ~ bbs(movingSpeed) + bbs(compassDirection) +
                                        bspatial(latitude, longitude) + 
                                        bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) + 
                                        bols(nrStatus, contrasts.arg = "contr.sum") +
                                        bbs(lte_rssi) +  bbs(lte_rsrp) + bbs(lte_rsrq) +
                                        bbs(nr_ssRsrp) + bbs(nr_ssSinr), 
                                    data = train_df, 
                                    na.action = na.pass,
                                    control = boost_control(mstop = mstop_init))
                
                cv_gam_train = cvrisk(gam_train, folds = cv(weights = model.weights(gam_train), type = "kfold", B = inner_folds))
                gam_train[mstop(cv_gam_train)]
                gam_preds = gam_train$predict(test_df)
                
                ## Black box
                bb_train = gbm(Throughput ~ movingSpeed + compassDirection +
                                latitude + longitude +
                                Throughput_lag1 + Throughput_lag2 + Throughput_lag3 + 
                                nrStatus +
                                lte_rssi +  lte_rsrp + lte_rsrq +
                                nr_ssRsrp + nr_ssSinr, 
                            data = train_df, 
                            distribution = "gaussian",
                            n.trees = mstop_init,
                            cv.folds = inner_folds,
                            n.cores = 1)
                
                #Tuning the number of boosting iterations
                bb_iter = gbm.perf(bb_train, plot.it = F, method = "cv")
                bb_preds = predict(bb_train, newdata = test_df, n.trees = bb_iter)
                
                ## Complete cases (CC) GAM
                cc_gam_train = gamboost(Throughput ~ bbs(movingSpeed) + bbs(compassDirection) + #M
                                            bspatial(latitude, longitude) + #L
                                            bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) + 
                                            bbs(lte_rssi) +  bbs(lte_rsrp) + bbs(lte_rsrq) +
                                            bbs(nr_ssRsrp) + bbs(nr_ssSinr), 
                                        data = train_df, 
                                        na.action = na.omit,
                                        control = boost_control(mstop = mstop_init))
                
                cv_cc_gam = cvrisk(cc_gam_train, folds = cv(weights = model.weights(cc_gam_train), type = "kfold", B = inner_folds))
                cc_gam_train[mstop(cv_cc_gam)]
                cc_gam_preds = cc_gam_train$predict(test_df)
                
                ## Complete predictors(CP) GAM
                non_na_predictors = names(which(apply(train_df, 2, function(c) sum(is.na(c)) ) == 0))
                cp_mod_predictors = setdiff(non_na_predictors, c("run_num", "seq_num", "nrStatus"))
                
                cp_gam_train = gamboost(Throughput ~ .,
                                        baselearner = "bbs",
                                        data = train_df[, colnames(train_df) %in% cp_mod_predictors], 
                                        na.action = na.omit,
                                        control = boost_control(mstop = mstop_init))
                
                cv_cp_gam = cvrisk(cp_gam_train, folds = cv(weights = model.weights(cp_gam_train), type = "kfold", B = inner_folds))
                cp_gam_train[mstop(cv_cp_gam)]
                cp_gam_preds = cp_gam_train$predict(test_df)
                
                ## Imputed predictor (IMP) GAM
                imp_train_df = train_df
                imp_test_df = test_df
                
                for(i in 1:ncol(imp_train_df)){
                    if(is.numeric(imp_train_df[,i])){
                        var_mean = mean(imp_train_df[,i], na.rm = T)
                    } else if (is.factor(imp_train_df[,i])){
                        var_mean = names(which.max(table(imp_train_df[,i])))
                    }
                    imp_train_df[is.na(imp_train_df[,i]), i] = var_mean
                    imp_test_df[is.na(imp_test_df[,i]), i] = var_mean
                }
                
                imp_gam_train = gamboost(Throughput ~ bbs(movingSpeed) + bbs(compassDirection) + #M
                                            bspatial(latitude, longitude) + #L
                                            bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) + 
                                            bbs(lte_rssi) +  bbs(lte_rsrp) + bbs(lte_rsrq) +
                                            bbs(nr_ssRsrp) + bbs(nr_ssSinr), 
                                        data = imp_train_df, 
                                        na.action = na.omit,
                                        control = boost_control(mstop = mstop_init))
                
                cv_imp_gam = cvrisk(imp_gam_train, folds = cv(weights = model.weights(imp_gam_train), type = "kfold", B = inner_folds))
                imp_gam_train[mstop(cv_imp_gam)]
                imp_gam_preds = imp_gam_train$predict(imp_test_df)
                
                ## LM fit with BAMMO (LM)
                lm_train = gamboost(Throughput ~ bols(movingSpeed) + bols(compassDirection) +
                                            bols(latitude, longitude) + 
                                            bols(Throughput_lag1) + bols(Throughput_lag2) + bols(Throughput_lag3) + 
                                            bols(nrStatus, contrasts.arg = "contr.sum") +
                                            bols(lte_rssi) +  bols(lte_rsrp) + bols(lte_rsrq) +
                                            bols(nr_ssRsrp) + bols(nr_ssSinr), 
                                        data = train_df, 
                                        na.action = na.pass,
                                        control = boost_control(mstop = mstop_init))
                
                cv_lm_train = cvrisk(lm_train, folds = cv(weights = model.weights(lm_train), type = "kfold", B = inner_folds))
                lm_train[mstop(cv_lm_train)]
                lm_preds = lm_train$predict(test_df)
                
                ## Baseline AR(1) model
                ar = lm(Throughput ~ Throughput_lag1, 
                        data = train_df)
                
                ar_preds = predict(ar, test_df)
                
                ## 3 lags model
                lag3_train = gamboost(Throughput ~  bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3), 
                                        data = train_df, 
                                        na.action = na.omit,
                                        control = boost_control(mstop = mstop_init))
                
                cv_lag3_train = cvrisk(lag3_train, folds = cv(weights = model.weights(lag3_train), type = "kfold", B = inner_folds))
                lag3_train[mstop(cv_lag3_train)]
                lag3_preds = lag3_train$predict(test_df)
                
                #MICE imputation
                mice_train <- mice(train_df[, !(colnames(train_df) %in% rm_vars)], printFlag = F, m = mice_imps)
                mice_test <- mice(test_df[, !(colnames(train_df) %in% rm_vars)], printFlag = F, m = mice_imps)
                
                mice_pred_mat <- matrix(NA, nrow = nrow(test_df), ncol = mice_imps)
                
                for(m in 1:mice_imps){
                    mice_train_m <- data.frame(Throughput = train_df$Throughput, complete(mice_train, m))

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
                
                mice_gam_preds <- apply(mice_pred_mat, 1, mean)


                #missForest imputation via missRanger
                lumos_missRanger_train <- missRanger(train_df[, !(colnames(train_df) %in% rm_vars)], num.trees = 50, data_only = T, verbose = F)
            
                missRanger_gam_train <- gamboost(train_df$Throughput ~ bbs(movingSpeed) + bbs(compassDirection) + #M
                                                    bspatial(latitude, longitude) + #L
                                                    bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) + 
                                                    bbs(lte_rssi) +  bbs(lte_rsrp) + bbs(lte_rsrq) +
                                                    bbs(nr_ssRsrp) + bbs(nr_ssSinr),
                                                data = lumos_missRanger_train,
                                                na.action = na.omit,
                                                control = boost_control(mstop = mstop_init))
                
                cv_missRanger_gam <- cvrisk(missRanger_gam_train, folds = cv(weights = model.weights(missRanger_gam_train), type = "kfold", B = inner_folds))
                missRanger_gam_train[mstop(cv_missRanger_gam)]
                
                lumos_missRanger_test <- missRanger(test_df[, !(colnames(test_df) %in% rm_vars)], num.trees = 50, data_only = T, verbose = F)
                missRanger_gam_preds <- missRanger_gam_train$predict(lumos_missRanger_test)

                # Make predictions on test data and compute MSE and MAE
                mses = rep(NA, 10)
                maes = rep(NA, 10)
                
                ## BAMMO
                mses[1] = mean((test_df$Throughput - gam_preds)^2)
                maes[1] = mean(abs(test_df$Throughput - gam_preds))
                
                ## Black box
                mses[2] = mean((test_df$Throughput - bb_preds)^2)
                maes[2] = mean(abs(test_df$Throughput - bb_preds))
                
                ## CC GAM
                mses[3] = mean((test_df$Throughput - cc_gam_preds)^2)
                maes[3] = mean(abs(test_df$Throughput - cc_gam_preds))
                
                ## CP GAM
                mses[4] = mean((test_df$Throughput - cp_gam_preds)^2)
                maes[4] = mean(abs(test_df$Throughput - cp_gam_preds))
                
                ## IMP GAM
                mses[5] = mean((test_df$Throughput - imp_gam_preds)^2)
                maes[5] = mean(abs(test_df$Throughput - imp_gam_preds))
                
                ## LM
                mses[6] = mean((test_df$Throughput - lm_preds)^2)
                maes[6] = mean(abs(test_df$Throughput - lm_preds))
                
                ## AR(1)
                mses[7] = mean((test_df$Throughput - ar_preds)^2)
                maes[7] = mean(abs(test_df$Throughput - ar_preds))

                ## lag3 model
                mses[8] = mean((test_df$Throughput - lag3_preds)^2)
                maes[8] = mean(abs(test_df$Throughput - lag3_preds))
                
                #MICE GAM
                mses[9] = mean((test_df$Throughput - mice_gam_preds)^2)
                maes[9] = mean(abs(test_df$Throughput - mice_gam_preds))

                #missRanger GAM
                mses[10] <- mean((test_df$Throughput - missRanger_gam_preds)^2)
                maes[10] <- mean(abs(test_df$Throughput - missRanger_gam_preds))
        
                FALSE
            },
            error = function(e){
                if(try_iter > 4) stop("too many tries!")
                print(e)
                TRUE
            })

            if(tryagain == FALSE){
                break
            }
            
            try_iter = try_iter + 1
        }

        return(
            list(
                mses = mses,
                maes = maes
            )
        )
    }
    
    #Run CV
    loop_out = mclapply(1:reps, sim_loop, mc.cores = ncores)
    
    #Collate output
    mses = t( sapply(loop_out, function(l) l$mses) )
    maes = t( sapply(loop_out, function(l) l$maes) )
    
    mods = c("GAM", "BB", "CC GAM", "CP GAM", "IMP GAM",  "LM", "AR1", "3 Lag", "MICE", "missRanger")
    
    colnames(mses) = mods
    colnames(maes) = mods
    
    #Create nice summary
    cv_summary = matrix(NA, nrow = length(mods), ncol = 4)
    colnames(cv_summary) = c("RMSE mean", "RMSE se", "MAE mean", "MAE se")
    rownames(cv_summary) = mods
    
    cv_summary[,1] = apply(sqrt(mses), 2, mean)
    cv_summary[,2] = apply(sqrt(mses), 2, function(c) sd(c)/sqrt( length(c) ) )
    cv_summary[,3] = apply(maes, 2, mean)
    cv_summary[,4] = apply(maes, 2, function(c) sd(c)/sqrt( length(c) ) )
    
    settings <- list(reps = reps, test_prop = test_prop, inner_folds = inner_folds, mstop_init = mstop_init, 
                     train_na_vars = train_na_vars, test_na_vars = test_na_vars, train_na_prop = train_na_prop, test_na_prop = test_na_prop, 
                     mice_imps = mice_imps)
    
    return(
        list(MAE = maes,
             MSE = mses,
             cv_summary = cv_summary,
             settings = settings)
    )
    
}

#### Run on Lumos data ####
cv_out_lumos = lumos_cv(df = lumos, reps = 100, test_prop = 1/10, inner_folds = 10, mstop_init = 100,
                        train_na_vars = NULL, test_na_vars = NULL, 
                        train_na_prop = 0, test_na_prop = 0, mice_imps = 5)

save(list = c("cv_out_lumos"), file = "lumos_cv_comparison.rds")

#### Run L,M,C scenarios ####

# Define groups of vars
L = c("latitude", "longitude")
M = c("movingSpeed", "compassDirection")
C = c("lte_rssi", "lte_rsrp", "lte_rsrq", "nr_ssRsrp", "nr_ssSinr")

# Iterate through combinations of (train_na_vars, test_na_vars)
train_na_var_list = list(NULL, L, M)
test_na_var_list = list(NULL, L, M)

cv_out_list = vector("list", length = length(train_na_var_list))

for(i in 1:length(train_na_var_list)){
    cv_out_list[[i]] = vector("list", length = length(test_na_var_list))
    for(j in 1:length(test_na_var_list)){
        cv_out_list[[i]][[j]] = lumos_cv(df = lumos, reps = 100, test_prop = 1/10, inner_folds = 10, mstop_init = 100,
                                         train_na_vars = train_na_var_list[[i]], test_na_vars = test_na_var_list[[j]],
                                         train_na_prop = 0.25, test_na_prop = 0.25, mice_imps = 5)
    }
}

# Save output
save(list = c("cv_out_list", "train_na_var_list", "test_na_var_list"), file = "lumos_cv_comparison_na_scenarios.rds")

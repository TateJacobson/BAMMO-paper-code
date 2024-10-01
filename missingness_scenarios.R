library(mboost)
library(gbm)
library(mice)
library(missRanger)
library(parallel)
library(missMethods)

ppn <- as.numeric(Sys.getenv("SLURM_NTASKS"))
ncores <- ifelse(is.na(ppn), 1, ppn)

#Run script to load and clean data
source('lumos_cleaning.R')

# Function for MCAR, MAR, MNAR prediction comparison  #####
lumos_missing_sims <- function(df, reps = 100, test_prop = 1/10, inner_folds = 10, mstop_init = 100, 
                               type = c("MCAR", "MAR", "MNAR"), miss_prop = 0.25, cols_mis = NULL, cols_ctrl = NULL, or = NULL,
                               grouped = T, miss_only = T, mice_imps = 5){
    
    type <- match.arg(type)
    
    # Variables to exclude from imputation
    rm_vars <- c("run_num", "seq_num", "lte_rssnr", "nr_ssRsrq", "Throughput", "nrStatus")
    
    sim_loop <- function(r){

        tryagain = FALSE
        try_iter = 1

        repeat {
            tryagain = tryCatch({
                # generate NAs using mechanism specified by type
                if(type == "MCAR"){
                    if(grouped){
                        dfn <- delete_MCAR(df, p = miss_prop, cols_mis = cols_mis[1])
                        dfn[ is.na( dfn[, cols_mis[1]] ) ,cols_mis[-1]] <- NA
                    } else {
                        dfn <- delete_MCAR(df, p = miss_prop, cols_mis = cols_mis)
                    }
                } else if (type == "MAR"){
                    if(is.null(cols_mis) || is.null(cols_ctrl) || is.null(or)) stop("missing necessary args for MAR")
                    if(grouped){
                        dfn <- delete_MAR_1_to_x(df, p = miss_prop, cols_mis = cols_mis[1], cols_ctrl = cols_ctrl[1], x = or)
                        dfn[ is.na( dfn[, cols_mis[1]] ) ,cols_mis[-1]] <- NA
                    } else {
                        dfn <- delete_MAR_1_to_x(df, p = miss_prop, cols_mis = cols_mis, cols_ctrl = cols_ctrl, x = or)
                    }
                    
                } else {
                    if(is.null(cols_mis) || is.null(or)) stop("missing necessary args for MAR")
                    if(grouped){
                        dfn <- delete_MNAR_1_to_x(df, p = miss_prop, cols_mis = cols_mis[1], x = or)
                        dfn[ is.na( dfn[,cols_mis[1]] ) ,cols_mis[-1]] <- NA
                    } else {
                        dfn <- delete_MNAR_1_to_x(df, p = miss_prop, cols_mis = cols_mis, x = or)
                    }
                }
                
                # pull training and test data
                n <- nrow(dfn)
                test_inds <- sample(n, n*test_prop)
                
                test_df <- dfn[test_inds,]
                train_df <- dfn[-test_inds,]
                
                # Fit models to training data and tune mstop
                ## BAMMO
                gam_train = gamboost(Throughput ~ bbs(movingSpeed) + bbs(compassDirection) +
                                        bspatial(latitude, longitude) +
                                        bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) +
                                        bbs(lte_rssi) +  bbs(lte_rsrp) + bbs(lte_rsrq) +
                                        bbs(nr_ssRsrp) + bbs(nr_ssSinr),
                                    data = train_df,
                                    na.action = na.pass,
                                    control = boost_control(mstop = mstop_init))

                cv_gam_train = cvrisk(gam_train, folds = cv(weights = model.weights(gam_train), type = "kfold", B = inner_folds))
                gam_train[mstop(cv_gam_train)]
                gam_preds = gam_train$predict(test_df)

                if(!miss_only){
                    ## Black box
                    bb_train = gbm(Throughput ~ movingSpeed + compassDirection +
                                    latitude + longitude +
                                    Throughput_lag1 + Throughput_lag2 + Throughput_lag3 +
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
                }
                
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
                cp_mod_predictors = setdiff(non_na_predictors, c("run_num", "seq_num", "nrStatus", "lte_rssnr", "nr_ssRsrq"))

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
                
                if(!miss_only){
                    ## LM fit with BAMMO (LM)
                    lm_train = gamboost(Throughput ~ bols(movingSpeed) + bols(compassDirection) +
                                            bols(latitude, longitude) + 
                                            bols(Throughput_lag1) + bols(Throughput_lag2) + bols(Throughput_lag3) + 
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
                }
                
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
                
                if(!miss_only){
                    ## Black box
                    mses[2] = mean((test_df$Throughput - bb_preds)^2)
                    maes[2] = mean(abs(test_df$Throughput - bb_preds))
                }
                
                ## CC GAM
                mses[3] = mean((test_df$Throughput - cc_gam_preds)^2)
                maes[3] = mean(abs(test_df$Throughput - cc_gam_preds))
                
                ## CP GAM
                mses[4] = mean((test_df$Throughput - cp_gam_preds)^2)
                maes[4] = mean(abs(test_df$Throughput - cp_gam_preds))
                
                ## IMP GAM
                mses[5] = mean((test_df$Throughput - imp_gam_preds)^2)
                maes[5] = mean(abs(test_df$Throughput - imp_gam_preds))
                
                if(!miss_only){
                    ## LM
                    mses[6] = mean((test_df$Throughput - lm_preds)^2)
                    maes[6] = mean(abs(test_df$Throughput - lm_preds))
                    
                    ## AR(1)
                    mses[7] = mean((test_df$Throughput - ar_preds)^2)
                    maes[7] = mean(abs(test_df$Throughput - ar_preds))
                    
                    ## lag3 model
                    mses[8] = mean((test_df$Throughput - lag3_preds)^2)
                    maes[8] = mean(abs(test_df$Throughput - lag3_preds))
                }
                
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
                     type = type, miss_prop = miss_prop, cols_mis = cols_mis, cols_ctrl = cols_ctrl, or = or, 
                     mice_imps = mice_imps)
    
    return(
        list(MAE = maes,
             MSE = mses,
             cv_summary = cv_summary,
             settings = settings)
    )
}

# Run MCAR, MAR, MNAR scenarios #####

#get only complete cases
lumos_complete <- lumos[complete.cases(lumos),]

mcar_out <- lumos_missing_sims(lumos_complete, reps = 100, test_prop = 1/10, inner_folds = 10, mstop_init = 100, 
                               type = "MCAR", miss_prop = 0.25, 
                               cols_mis = c("nr_ssRsrp","nr_ssSinr"),
                               grouped = T, miss_only = T, mice_imps = 5)

save(list = c("mcar_out"), file = "lumos_MCAR_sims.rds")
 
mar_out <- lumos_missing_sims(lumos_complete, reps = 100, test_prop = 1/10, inner_folds = 10, mstop_init = 100, 
                              type = "MAR", miss_prop = 0.25, 
                              cols_mis= c("nr_ssRsrp","nr_ssSinr"), 
                              cols_ctrl = c("movingSpeed", "movingSpeed"), or = 9,
                              grouped = T, miss_only = T, mice_imps = 5)
 
save(list = c("mar_out"), file = "lumos_MAR_sims.rds")
 
mnar_out <- lumos_missing_sims(lumos_complete, reps = 100, test_prop = 1/10, inner_folds = 10, mstop_init = 100,
                               type = "MNAR", miss_prop = 0.25,
                               cols_mis= c("nr_ssRsrp","nr_ssSinr"), or = 1/9,
                               grouped = T, miss_only = T, mice_imps = 5)
 
save(list = c("mnar_out"), file = "lumos_MNAR_sims.rds")
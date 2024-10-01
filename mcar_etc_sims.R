library(mboost)
library(gbm)
library(mice)
library(missRanger)
library(parallel)
library(missMethods)
library(MASS)

ppn <- as.numeric(Sys.getenv("SLURM_NTASKS"))
ncores <- ifelse(is.na(ppn), 1, ppn)

# Function for MCAR, MAR, MNAR prediction comparison  #####
lumos_missing_sims <- function(ntrain, ntest, p, rho, beta, 
                               reps = 100, inner_folds = 10, mstop_init = 100, 
                               type = c("MCAR", "MAR", "MNAR"), miss_prop = 0.25, cols_mis = NULL, cols_ctrl = NULL, or = NULL,
                               grouped = T, mice_imps = 5){
    
    type <- match.arg(type)
    
    n <- ntrain + ntest
    
    sim_loop <- function(r){
        
        tryagain = FALSE
        try_iter = 1
        
        repeat {
            tryagain = tryCatch({

                #generate data
                SigmaX <- outer(1:p, 1:p, function(i, j) rho^abs(i-j))
                x <- mvrnorm(n, rep(0,p), SigmaX)
                
                y <- x%*%beta + rnorm(n, 0, 1)
                
                df <- data.frame(y,x)
                
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
                    
                    if(is.null(cols_mis) || is.null(or)) stop("missing necessary args for MNAR")
                    
                    if(grouped){
                        dfn <- delete_MNAR_1_to_x(df, p = miss_prop, cols_mis = cols_mis[1], x = or)
                        dfn[ is.na( dfn[,cols_mis[1]] ) ,cols_mis[-1]] <- NA
                    } else {
                        dfn <- delete_MNAR_1_to_x(df, p = miss_prop, cols_mis = cols_mis, x = or)
                    }
                    
                }
                
                # split into training and test data
                test_inds <- sample(n, ntest)
                
                test_df <- dfn[test_inds,]
                train_df <- dfn[-test_inds,]
                
                # Fit models to training data and tune mstop
                ## BAMMO
                gam_train = gamboost(y ~ .,
                                     baselearner = "bbs",
                                     data = train_df,
                                     na.action = na.pass,
                                     control = boost_control(mstop = mstop_init))
                
                cv_gam_train = cvrisk(gam_train, folds = cv(weights = model.weights(gam_train), type = "kfold", B = inner_folds))
                gam_train[mstop(cv_gam_train)]
                gam_preds = gam_train$predict(test_df)
                
                
                ## Complete cases (CC) GAM
                cc_gam_train = gamboost(y ~ .,
                                        baselearner = "bbs",
                                        data = train_df,
                                        na.action = na.omit,
                                        control = boost_control(mstop = mstop_init))
                
                cv_cc_gam = cvrisk(cc_gam_train, folds = cv(weights = model.weights(cc_gam_train), type = "kfold", B = inner_folds))
                cc_gam_train[mstop(cv_cc_gam)]
                cc_gam_preds = cc_gam_train$predict(test_df)
                
                ## Complete predictors(CP) GAM
                non_na_predictors <- which( apply(train_df, 2, function(c) sum(is.na(c)) ) == 0)
                
                cp_gam_train <- gamboost(y ~ .,
                                         baselearner = "bbs",
                                         data = train_df[, non_na_predictors],
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
                
                imp_gam_train <- gamboost(y ~ .,
                                          baselearner = "bbs",
                                          data = imp_train_df,
                                          na.action = na.omit,
                                          control = boost_control(mstop = mstop_init))
                
                cv_imp_gam = cvrisk(imp_gam_train, folds = cv(weights = model.weights(imp_gam_train), type = "kfold", B = inner_folds))
                imp_gam_train[mstop(cv_imp_gam)]
                imp_gam_preds = imp_gam_train$predict(imp_test_df)
                
                #MICE imputation
                mice_train <- mice(train_df[, -1], printFlag = F, m = mice_imps)
                mice_test <- mice(test_df[, -1], printFlag = F, m = mice_imps)
                
                mice_pred_mat <- matrix(NA, nrow = ntest, ncol = mice_imps)
                
                for(m in 1:mice_imps){
                    mice_train_m <- data.frame(y = train_df$y, complete(mice_train, m))
                    
                    mice_gam_train <- gamboost(y ~ .,
                                               baselearner = "bbs",
                                               data = mice_train_m,
                                               na.action = na.omit,
                                               control = boost_control(mstop = mstop_init))
                    
                    cv_mice_gam <- cvrisk(mice_gam_train, folds = cv(weights = model.weights(mice_gam_train), type = "kfold", B = inner_folds))
                    mice_gam_train[mstop(cv_mice_gam)]
                    
                    mice_pred_mat[, m] <- mice_gam_train$predict(complete(mice_test, m)) 
                }
                
                mice_gam_preds <- apply(mice_pred_mat, 1, mean)
                
                #missForest imputation via missRanger
                missRanger_train <- missRanger(train_df[, -1], num.trees = 50, data_only = T, verbose = F)
                missRanger_test <- missRanger(test_df[, -1], num.trees = 50, data_only = T, verbose = F)
                
                missRanger_train <- data.frame(y = train_df$y, missRanger_train)
                
                missRanger_gam_train <- gamboost(y ~ .,
                                                 baselearner = "bbs",
                                                 data = missRanger_train,
                                                 na.action = na.omit,
                                                 control = boost_control(mstop = mstop_init))
                
                cv_missRanger_gam <- cvrisk(missRanger_gam_train, folds = cv(weights = model.weights(missRanger_gam_train), type = "kfold", B = inner_folds))
                missRanger_gam_train[mstop(cv_missRanger_gam)]
                
                missRanger_gam_preds <- missRanger_gam_train$predict(missRanger_test)
                
                # Make predictions on test data and compute MSE and MAE
                mses = rep(NA, 6)
                maes = rep(NA, 6)
                
                ## Our GAM
                mses[1] = mean((test_df$y - gam_preds)^2)
                maes[1] = mean(abs(test_df$y - gam_preds))
                
                ## CC GAM
                mses[2] = mean((test_df$y - cc_gam_preds)^2)
                maes[2] = mean(abs(test_df$y - cc_gam_preds))
                
                ## CP GAM
                mses[3] = mean((test_df$y - cp_gam_preds)^2)
                maes[3] = mean(abs(test_df$y - cp_gam_preds))
                
                ## IMP GAM
                mses[4] = mean((test_df$y - imp_gam_preds)^2)
                maes[4] = mean(abs(test_df$y - imp_gam_preds))
                
                #MICE GAM
                mses[5] = mean((test_df$y - mice_gam_preds)^2)
                maes[5] = mean(abs(test_df$y - mice_gam_preds))
                
                #missRanger GAM
                mses[6] <- mean((test_df$y - missRanger_gam_preds)^2)
                maes[6] <- mean(abs(test_df$y - missRanger_gam_preds))
                
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
    loop_out <- mclapply(1:reps, sim_loop, mc.cores = ncores)
    
    #Collate output
    mses <- t( sapply(loop_out, function(l) l$mses) )
    maes <- t( sapply(loop_out, function(l) l$maes) )
    
    mods <- c("GAM", "CC GAM", "CP GAM", "IMP GAM", "MICE", "missRanger")
    
    colnames(mses) <- mods
    colnames(maes) <- mods
    
    #Create nice summary
    cv_summary <- matrix(NA, nrow = length(mods), ncol = 4)
    colnames(cv_summary) <- c("RMSE mean", "RMSE se", "MAE mean", "MAE se")
    rownames(cv_summary) <- mods
    
    cv_summary[,1] <- apply(sqrt(mses), 2, mean)
    cv_summary[,2] <- apply(sqrt(mses), 2, function(c) sd(c)/sqrt( length(c) ) )
    cv_summary[,3] <- apply(maes, 2, mean)
    cv_summary[,4] <- apply(maes, 2, function(c) sd(c)/sqrt( length(c) ) )
    
    settings <- list(ntrain = ntrain, ntest = ntest, p = p, rho = rho, beta = beta, 
                     reps = reps, inner_folds = inner_folds, mstop_init = mstop_init, 
                     type = type, miss_prop = miss_prop, cols_mis = cols_mis, cols_ctrl = cols_ctrl, or = or,
                     grouped = grouped)
    
    return(
        list(MAE = maes,
             MSE = mses,
             cv_summary = cv_summary,
             settings = settings)
    )
}

# Run MCAR, MAR, MNAR scenarios #####
ntrain <- 100 
ntest <- 100
p <- 10
rho <- 0.5
beta <-c(1, 0, 1/2, -2, rep(0, p - 4))

mcar_out <- lumos_missing_sims(ntrain = ntrain, ntest = ntest, p = p, rho = rho, beta = beta,
                               reps = 100, inner_folds = 10, mstop_init = 100, 
                               type = "MCAR", miss_prop = 0.25, 
                               cols_mis = c(2,3),
                               grouped = F, mice_imps = 5)

save(list = c("mcar_out"), file = "MCAR_etc_sims.rds")

mar_out <- lumos_missing_sims(ntrain = ntrain, ntest = ntest, p = p, rho = rho, beta = beta,
                              reps = 100, inner_folds = 10, mstop_init = 100, 
                              type = "MAR", miss_prop = 0.25, 
                              cols_mis = c(2, 3),
                              cols_ctrl = c(4, 4), or = 9,
                              grouped = F, mice_imps = 5)

save(list = c("mcar_out", "mar_out"), file = "MCAR_etc_sims.rds")

mnar_out <- lumos_missing_sims(ntrain = ntrain, ntest = ntest, p = p, rho = rho, beta = beta,
                               reps = 100, inner_folds = 10, mstop_init = 100, 
                               type = "MNAR", miss_prop = 0.25, 
                               cols_mis = c(2,3), or = 1/9,
                               grouped = F, mice_imps = 5)

save(list = c("mcar_out", "mar_out", "mnar_out"), file = "MCAR_etc_sims.rds")

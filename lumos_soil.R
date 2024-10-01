library(mboost)
library(tidyverse)

#Run script to load and clean data
source("lumos_cleaning.R")

#Fit initial BAMMO model
gammo_form = Throughput ~ bbs(movingSpeed) + bbs(compassDirection) + 
    bspatial(latitude, longitude) + 
    bbs(Throughput_lag1) + bbs(Throughput_lag2) + bbs(Throughput_lag3) + 
    bols(nrStatus, contrasts.arg = "contr.sum") +
    bbs(lte_rssi) +  bbs(lte_rsrp) + bbs(lte_rsrq) +
    bbs(nr_ssRsrp) + bbs(nr_ssSinr)  

gammo_control = boost_control(mstop = 100, nu = 0.1)

gammo = gamboost(gammo_form,
                 data = lumos, 
                 na.action = na.pass,
                 control = gammo_control)

##### Weighting function #####

# Function to compute dfs for all selected base learners
compute_df =  function(object){
    nu = object$control$nu
    xs = object$xselect()
    vars = sort(unique(xs))
    bl_names = names(object$baselearner)
    
    w = object$`(weights)`
    bl_dfs = sapply(vars, function(v)  object$baselearner[[v]]$dpp(w)$df()[[1]])
    bl_dfTrS = sapply(vars, function(v)  object$baselearner[[v]]$dpp(w)$df()[[3]])
    
    names(bl_dfs) = bl_names[vars]
    return(list(dfs = bl_dfs, dfTrS = bl_dfTrS ))
}

dfs <- compute_df(gammo)

#Function to compute ARM weights
compute_weights = function(object, dat, psi = 1, 
                           L = 10, test_prop = 0.1, 
                           control = boost_control(mstop = 100, nu = 0.1)){
    
    M = object$mstop()
    n =  nrow(dat)
    
    #Get sequence of selected variables
    vars = names(object$baselearner)
    xs = object$xselect()
    vars_selected = vars[xs]
    
    arm_wts = matrix(NA, ncol = M, nrow = L)
    
    for(l in 1:L){
        #Split data
        test_inds = sample(1:n, n*test_prop)
        test_dat = dat[test_inds,]
        train_dat = dat[-test_inds,]
        
        loglik_test = rep(NA, M)
        dfs = rep(NA, M)
        
        #Iterate through models A_1,...,A_M
        for(m in 1:M){
            #Set variables and formula to use in mth model
            vars_m = unique(vars_selected[1:m])
            if(m > 1 & (length(vars_m) == length( unique( vars_selected[1:(m-1)] ) ) )  ){
                dfs[m] = dfs[m-1]
                loglik_test[m] = loglik_test[m-1]
            } else {
                form_m = reformulate(termlabels = vars_m, response = "Throughput")
                gam_m = gamboost(form_m, data = train_dat, na.action = na.pass, control = control)
                
                dfs[m] = sum( compute_df(gam_m)$dfTrS )
                s2_m = sum((gam_m$response - gam_m$fitted())^2)/(length(gam_m$response) - dfs[m] - 1)
                
                if(nrow(test_dat) > 1e4){ #have to chunk test data for prediction if test data set has too many observations
                    n_test_chunks = nrow(test_dat)/1e4
                    yhat_m = c()
                    for(tc in 1:n_test_chunks){
                        start_ind = 1e4*(tc-1) + 1
                        if(tc == n_test_chunks){
                            end_ind = nrow(test_dat)
                        } else {
                            end_ind = 1e4*tc
                        }
                        test_chunk = test_dat[start_ind:end_ind,]  
                        yhat_m = c(yhat_m, gam_m$predict(test_chunk))
                    }
                } else {
                    yhat_m = gam_m$predict(test_dat)
                }
                
                loglik_test[m] = -sum((test_dat$Throughput - yhat_m)^2)/(2*s2_m) - n/2*log(s2_m)
            }
        }
        
        C_prior = dfs*log( exp(1)*max(dfs)/dfs ) + 2*log(dfs + 2)
        
        for(m in 1:M){
            loglik_diffs_m = loglik_test - loglik_test[m]
            prior_diffs_m = C_prior - C_prior[m]
            arm_wts[l, m] = 1/(sum( exp(loglik_diffs_m - psi*prior_diffs_m) ))
        }
    }
    
    wts = apply(arm_wts, 2, mean)
    
    
    return(wts)
}

gammo_arm_wts <- compute_weights(gammo, dat = lumos, L = 100, test_prop = 0.5)

##### Compute SOIL #####

compute_soil = function(object, weights){
    M = object$mstop()
    
    #Get sequence of selected variables
    vars = names(object$baselearner)
    xs = object$xselect()
    vars_selected = vars[xs]
    
    first_selection = rep(NA, length(vars))
    soil = rep(0, length(vars))
    for(i in 1:length(vars)){
        #Check if vars[i] is in each model
        in_Am = rep(0, M)
        if(vars[i] %in% vars_selected){
            first_selection[i] = which(vars_selected == vars[i])[1]
            in_Am[(first_selection[i]):M] = 1
            soil[i] = sum(in_Am*weights)
        }  
    }
    
    names(soil) = vars
    return(soil)
}

soil_arm = compute_soil(gammo, gammo_arm_wts)

##### Plot SOIL output #####
var_names =  c("moving speed", "compass direction", 
               "(latitude, longitude)", 
               "throughput (lag 1)", "throughput (lag 2)", "throughput (lag 3)", 
               "5G NR Status",
               "LTE RSSI", "LTE RSRP", "LTE RSRQ",
               "5G NR RSRP", "5G NR SSNR")

soil_df = data.frame(names = var_names, ARM = soil_arm)
soil_df = soil_df[order(soil_df$ARM, decreasing = T),]

soil_df %>%
    ggplot(mapping = aes(x = names, y = ARM)) +
    geom_col() + 
    scale_x_discrete(limits = rev(soil_df$names) ) + 
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    coord_flip() + 
    theme_bw() +
    theme(plot.margin = margin(10, 20, 10, 10, "points")) +
    labs(y = "SOIL Importance", x = "Variable")

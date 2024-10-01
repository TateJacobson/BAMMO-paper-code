library(ggplot2)
library(gridExtra)

# MSE and MAE side-by-side on original Lumos data ####
load("lumos_cv_comparison.rds")
og_out <- cv_out_lumos$cv_summary
method_names <- c("BAMMO", "BB", "CC", "CP", "IMP",  "LM", "AR1", "3 Lag", "MICE", "missForest")
og_out <- data.frame(method = method_names, og_out)

(og_out_miss <- og_out[c(1,3,4,5,9,10),])

(og_out_spec <- og_out[c(1,2,6,7,8),])

makeplots <- function(df, filename){
    shape <- 20
    
    rmse_plot <- ggplot(data = df) + 
        geom_pointrange(aes(x = method, y = RMSE.mean, ymax = RMSE.mean + RMSE.se, ymin = RMSE.mean - RMSE.se), shape = shape)+
        scale_x_discrete(limits = df$method) +
        labs(x = "Method", y = "RMSE")+
        theme_bw()
    
    mae_plot <- ggplot(data = df) + 
        geom_pointrange(aes(x = method, y = MAE.mean, ymax = MAE.mean + MAE.se, ymin = MAE.mean - MAE.se), shape = shape)+
        scale_x_discrete(limits = df$method) +
        labs(x = "Method", y = "MAE")+
        theme_bw()
    
    grid.arrange(grobs = list(rmse_plot, mae_plot), ncol = 2)
    
    pdf(filename, width = 8, height = 4)
    grid.arrange(grobs = list(rmse_plot, mae_plot), ncol = 2)
    dev.off()
}

makeplots(og_out_miss, "lumos_cv_miss.pdf")
makeplots(og_out_spec, "lumos_cv_spec.pdf")

# RMSE plots from L,M,C scenarios ####
load("lumos_cv_comparison_na_scenarios.rds")
shape <- 20
method_names = c("BAMMO", "BB", "CC", "CP", "IMP",  "LM", "AR1", "3 Lag", "MICE", "missForest")

make_na_scenario_plots <- function(type = c("miss", "spec"), metric = c("rmse", "mae")){
    
    plot_list = list()
    train_list = c("C", "C+L", "C+M")
    test_list = c("C", "C+L", "C+M")
    
    for(tr in 1:3){
        for(tst in 1:3){
            out = cv_out_list[[tr]][[tst]]$cv_summary
            out = data.frame(method = method_names, out)
            
            if(type == "miss"){
                out <- out[c(1,3,4,5,9,10),]
                
                rmse_y_lims <- c(172.5, 178)
                mae_y_lims <- c(106, 116)
                    
            } else if (type == "spec"){
                out <- out[c(1,2,6,7,8),]   
                
                rmse_y_lims <- c(172.5, 179.5)
                mae_y_lims <- c(106.5, 110.5)
            }
            
            if(metric == "rmse"){
                
                plot_list[[ 3*(tr-1) + tst ]] <- ggplot(data = out) + 
                    geom_pointrange(aes(x = method, y = RMSE.mean, ymax = RMSE.mean + RMSE.se, ymin = RMSE.mean - RMSE.se), shape = shape)+
                    scale_x_discrete(limits = out$method) +
                    scale_y_continuous(limits = rmse_y_lims) +
                    labs(x="", y="", title = paste0(train_list[[tr]], " train, ", test_list[[tst]]," test")) +
                    theme_bw() +
                    theme(plot.margin = margin(2, 2, 2, 2, "points"), 
                          plot.title = element_text(size = 10, hjust = 0.5),
                          axis.text = element_text(size = 7))
                
            } else if (metric == "mae"){
                
                plot_list[[ 3*(tr-1) + tst ]] <- ggplot(data = out) + 
                    geom_pointrange(aes(x = method, y = MAE.mean, ymax = MAE.mean + MAE.se, ymin = MAE.mean - MAE.se), shape = shape)+
                    scale_x_discrete(limits = out$method) +
                    scale_y_continuous(limits = mae_y_lims) +
                    labs(x="", y="", title = paste0(train_list[[tr]], " train, ", test_list[[tst]]," test")) +
                    theme_bw() +
                    theme(plot.margin = margin(2, 2, 2, 2, "points"), 
                          plot.title = element_text(size = 10, hjust = 0.5),
                          axis.text = element_text(size = 7))
                
            }
        }
    }
    
    filename = paste0(metric, "_na_prediction_comparison_", type, ".pdf")
    
    grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)
    
    pdf(filename, width = 10, height = 10)
    grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)
    dev.off()
}

make_na_scenario_plots("miss", "rmse")
make_na_scenario_plots("miss", "mae")
make_na_scenario_plots("spec", "rmse")
make_na_scenario_plots("spec", "mae")

# Plot results from MCAR, MAR, MNAR scenarios ####
make_miss_plots <- function(in_file, df_name, out_file, type = c("miss", "spec", "sims")){
    type <- match.arg(type)
    
    load(in_file)
    df <- get(df_name)
    
    og_out <- df$cv_summary
    
    if(type == "miss"){
        method_names <- c("BAMMO", "BB", "CC", "CP", "IMP",  "LM", "AR1", "3 Lag", "MICE", "missForest")
        og_out <- data.frame(method = method_names, og_out)
        og_out <- og_out[c(1,3,4,5,9,10),]
    } else if (type == "spec"){
        method_names <- c("BAMMO", "BB", "CC", "CP", "IMP",  "LM", "AR1", "3 Lag", "MICE", "missForest")
        og_out <- data.frame(method = method_names, og_out)
        og_out <- og_out[c(1,2,6,7,8),]
    } else if (type == "sims"){
        method_names <- c("BAMMO", "CC", "CP", "IMP", "MICE", "missForest")
        og_out <- data.frame(method = method_names, og_out)
    }
    
    makeplots(og_out, out_file)
}

# augmented Lumos data
make_miss_plots("lumos_MCAR_sims.rds", "mcar_out", out_file = "lumos_mcar.pdf", type = "miss")
make_miss_plots("lumos_MAR_sims.rds", "mar_out", out_file = "lumos_mar.pdf", type = "miss")
make_miss_plots("lumos_MNAR_sims.rds", "mnar_out", out_file = "lumos_mnar.pdf", type = "miss")

# simulated data 
make_miss_plots("MCAR_etc_sims.rds", "mcar_out", "mcar_sims.pdf", "sims")
make_miss_plots("MCAR_etc_sims.rds", "mar_out", "mar_sims.pdf", "sims")
make_miss_plots("MCAR_etc_sims.rds", "mnar_out", "mnar_sims.pdf", "sims")

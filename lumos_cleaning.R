##### Load and clean data #####

# Load Lumos5G data
lumos = read.csv("Lumos5G-v1.0/Lumos5G-v1.0.csv")

# Variable groupings:
## (L) Location based: latitude, longitude
## (M) Mobility based: movingSpeed, compassDirection
## (C) Connection based: past Throughput, nrStatus, lte_ variables, nr_ variables
lmc_cols = which(colnames(lumos) %in% c("Throughput", "latitude", "longitude", "movingSpeed", "compassDirection", "nrStatus",
                                        colnames(lumos)[startsWith(colnames(lumos), "lte_")],
                                        colnames(lumos)[startsWith(colnames(lumos), "nr_")]) )
index_cols = which(colnames(lumos) %in% c("run_num", "seq_num"))

# Remove unwanted variables
lumos_clean = lumos[,c(index_cols, lmc_cols)]
colnames(lumos_clean)

# Clean up variable types
lumos_clean$nrStatus = factor(lumos_clean$nrStatus)

# Function to create lagged Throughput columns for each run
lag_var = function(x_var, index_var, lags = 1){
    lag_list = list()
    for(l in 1:lags){
        lag_out = c()
        for(i in unique(index_var) ){
            xi = x_var[index_var == i]
            lag_out = c(lag_out, NA, xi[-length(xi)])
        }
        lag_list[[l]] = lag_out
        x_var = lag_out
    }
    names(lag_list) = paste0("Throughput_lag", 1:lags)
    data.frame(lag_list)
}

# Add lagged Throughput columns
lags = 3
lag_out = lag_var(lumos_clean$Throughput, lumos_clean$run_num, lags = lags)
lumos_clean2 = lumos_clean[lumos_clean$seq_num > lags,]
lumos_clean2 = cbind(lumos_clean2, lag_out[complete.cases(lag_out),])
head(lumos_clean2[lumos_clean2$run_num == 2, startsWith(colnames(lumos_clean2), "Throughput")])

rm(lumos)

lumos = lumos_clean2

# clean up
rm(lumos_clean, lumos_clean2, index_cols, lmc_cols, lag_out, lags, lag_var)

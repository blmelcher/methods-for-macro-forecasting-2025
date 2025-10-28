#in this file, we will test BVAR model with different priors, selected
#variables and lags to compute the RMSE in each case.

#results in this file will be used to compare different
# models and select the best one to use in main.R

#--------------------- libraries ---------------------#
library(ggplot2)
library(dplyr)
library(tidyr)
library(BVAR)
library(lubridate)
library(vars)  # needed for classical VAR benchmark

source("stationary.R")
# --------------------- set up ------------------------#
# set seed for reproducibility
set.seed(42)

# load data and data manipulation
df <- utils::read.csv("data/data_quarterly.csv")

df$date <- as.Date(paste0(df$date, "-01")) # format date
df <- df %>% filter(date <= as.Date("2025-07-01")) # until 01.10.2025

# inflate CPI to get inflation rate
df$inflation <- 100*(df$cpi - dplyr::lag(df$cpi, 1)) / dplyr::lag(df$cpi, 1) 
df <- df %>% filter(!is.na(inflation)) #remove first NA row

# gdp growth instead of gdp
df$gdp <- log(df$gdp)
# -------------------- apply log + growth transformations --------------------
# define the rate variables (do NOT log-transform these)
rate_variables <- c("inflation", "urilo", "srate", "srate_ge")

# safe log: add small offset if non-positive values present
#safe_log <- function(x, tiny = 1e-6) {
#  x_num <- as.numeric(x)
#  if (all(is.na(x_num))) return(x_num)
#  if (any(x_num <= 0, na.rm = TRUE)) {
#    offset <- abs(min(x_num, na.rm = TRUE)) + tiny
#    warning("Non-positive values found; adding offset = ", signif(offset, 6), " before log.")
#    x_num <- x_num + offset
#  }
#  return(log(x_num))
#}

# apply log to non-rate variables (excluding date)
# (var in names(df)) {
#  if (!(var %in% rate_variables) && var != "date") {
#    df[[var]] <- safe_log(df[[var]])
#  }
#}

# convert to growth rates (percent) for the same non-rate variables
for (var in names(df)) {
  if (!(var %in% rate_variables) && var != "date") {
    df[[var]] <-  (df[[var]] - dplyr::lag(df[[var]], 1))
  }
}

# remove initial rows with NA introduced by growth calculation if needed
# (subsequent code already removes rows with NA in gdp)

#check if stationary using functions from stationary.R
#for (var in setdiff(colnames(df), "date")) {  # excluding date column
#  cat("Checking stationarity for:", var, "\n")
#  ts_data <- ts(df[[var]])
#  stationary_ts <- make_stationary(ts_data, use_log = TRUE)   # returns aligned vector (same length)
#  df[[var]] <- as.numeric(stationary_ts)                     # assign directly, no extra NA
#}
df <- df %>% filter(!is.na(gdp))  # remove first row with NA

# ---------------------- lags ---------------------#

# possible lags: 1, 4 (one year), 8 (two years)
lags <- c(1,4)

# plot the forecast variables (save to files to avoid "figure margins too large")
fv <- c("gdp", "inflation", "wkfreuro")
dir.create("output/plots", recursive = TRUE, showWarnings = FALSE)
for (v in fv) {
  p <- ggplot(df, aes(x = date, y = .data[[v]])) +
    geom_line() +
    labs(title = v, x = "Date", y = v) +
    theme_bw()
  ggsave(
    filename = file.path("output/plots", paste0("plot_", v, ".png")),
    plot = p, width = 7, height = 4, dpi = 150
  )
}

# ------------------ selected variables ------------------#

#selected variables for BVAR model

selected_variables_0 <- c("gdp", "inflation", "wkfreuro")

selected_variables_1 <- c("gdp", "inflation", "wkfreuro", "consp", "consg", 
                          "ifix", "icnstr", "ime", "exc1", "imc1", "ltot", "uroff", "wage", 
                          "srate", "poilusd", "pcioecd", "vaabcde", "vaghji")

#from this initial set, we will also try a smaller set of variables

selected_variables_2 <- c("gdp", "inflation", "wkfreuro", "consp", "consg", "ifix", 
                          "exc1", "imc1", "ltot", "uroff", "wage", "srate", 
                          "poilusd", "pcioecd")

selected_variables_3 <- c("gdp", "inflation", "wkfreuro", "consp", "exc1", "ltot",
                          "wage", "srate")

#-------------------------------------------------------------------
# List of variable sets to test
variable_sets <- list(
  set_0 = selected_variables_0,
  set_1 = selected_variables_1,
  set_2 = selected_variables_2,
  set_3 = selected_variables_3
)

#--------------------- priors ---------------------#
# Define different prior configurations to test
mn <- bv_minnesota(
  lambda = bv_lambda(mode = 0.5, sd = 0.1, min = 0.001, max = 5),
  alpha  = bv_alpha(mode = 4),
)

soc <- bv_soc(mode = 1, sd = 0.5)   # shrink sum of AR coeffs to 1, with variance
priors <- bv_priors(
  hyper = c("lambda", "alpha", "psi"),
  mn = mn,
  soc = soc
)

# Diffuse (flat) prior: approximation to OLS using Minnesota with a very small lambda
diffuse_prior_bv <- bv_priors(
  hyper = c("lambda", "alpha"),
  mn = bv_minnesota(
    lambda = bv_lambda(mode = 1e-6, sd = 1e-6, min = 1e-12, max = 1),
    alpha  = bv_alpha(mode = 4)
  )
)

# Priors based on predefined components (mn and soc)
mn_prior_bv <- bv_priors(
  hyper = c("lambda", "alpha"),
  mn = mn
)

soc_prior_bv <- bv_priors(
  hyper = c("lambda", "alpha"),
  mn = mn,
  soc = soc
)

# List of priors ready to be passed to BVAR::bvar
prior_configs <- list(
  mn      = mn_prior_bv,
  soc     = soc_prior_bv
  # diffuse = diffuse_prior_bv
)

rolling_window_size <- 110  # e.g., 40 quarters (10 years)

#---------------------- end of setup ---------------------#
compute_bvar_rmse <- function(data, variables, lags, prior, window_size){
  # Rolling window BVAR RMSE, returns a data frame per variable
  horizon <- 1
  n_obs <- nrow(data)

  quantile_bands <- c("2.5%", "16%", "50%", "84%", "97.5%")
  forecast_variables <- c("gdp", "inflation", "wkfreuro")

  pred_q50 <- matrix(NA_real_, nrow = n_obs, ncol = length(variables),
                     dimnames = list(NULL, variables))
  
  for (i in seq(from = window_size + lags, to = n_obs - horizon)) {
    train_start <- i - window_size + 1
    train_end <- i
    
    y_train <- data[train_start:train_end, ]
    
    invisible(trained_model <- bvar(
      y_train %>% dplyr::select(all_of(variables)),
      lags = lags,
      n_draw = 10000,
      n_burn = 2500,
      n_thin = 1,
      priors = prior,
      verbose = FALSE
    ))
    
    prediction <- predict(trained_model, horizon = horizon, conf_bands = c(0.16, 0.025))
    
    pred_quants <- prediction$quants
    mat <- matrix(NA_real_, nrow = length(variables), ncol = length(quantile_bands),
                  dimnames = list(variables, quantile_bands))
    
    for (j in seq_along(variables)) {
      mat[j, ] <- pred_quants[, 1, j]
    }
    
    results <- as_tibble(mat, rownames = "variable") %>%
      rename(
        q025 = `2.5%`,
        q16  = `16%`,
        q50  = `50%`,
        q84  = `84%`,
        q975 = `97.5%`
      )
  
    pred_q50[i+1, ] <- results$q50
  }

  res <- data.frame(
    variable = forecast_variables,
    rmse = numeric(length(forecast_variables))
  )
  
  for (i in seq_along(forecast_variables)) {
    var <- forecast_variables[i]
    valid_indices <- which(!is.na(pred_q50[, var]))
    rmse <- sqrt(mean((pred_q50[valid_indices, var] - data[valid_indices, var])^2))
    res$rmse[i] <- rmse
  }
  
  return(res)
}


#–-------------------- BVAR/VAR RMSE computation ---------------------#

# data frame for results always add the new results at the end of the df, but also add columns with lags, variable set and prior used

results_df <- data.frame(
  variable = character(),
  rmse = numeric(),
  lags = integer(),
  variable_set = character(),
  prior = character(),
  stringsAsFactors = FALSE
)
for (lag in lags) {
  #print(paste("Testing lag:", lag))
  for (var_set_name in names(variable_sets)) {
    #cat("testing: ", var_set_name)
    selected_vars <- variable_sets[[var_set_name]]
    
    for (prior_name in names(prior_configs)) {
      cat("testing: ", lag, " ", var_set_name, " ",prior_name, "\n")
      prior_config <- prior_configs[[prior_name]]
      
      temp_results <- compute_bvar_rmse(
        data = df,
        variables = selected_vars,
        lags = lag,
        prior = prior_config,
        window_size = rolling_window_size
      )
      temp_results$lags <- lag
      temp_results$variable_set <- var_set_name
      temp_results$prior <- prior_name
      results_df <- rbind(results_df, temp_results)
      
    }
  }
}

results_df
# save results to the output folder
save(results_df, file = "output/results_df.RData")


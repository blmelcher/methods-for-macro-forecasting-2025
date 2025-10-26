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
df$inflation <- (df$cpi - dplyr::lag(df$cpi, 1)) / dplyr::lag(df$cpi, 1) 
df <- df %>% filter(!is.na(inflation)) #remove first NA row

# -------------------- apply log + growth transformations --------------------
# define the rate variables (do NOT log-transform these)
rate_variables <- c("inflation", "urilo", "srate", "srate_ge")

# safe log: add small offset if non-positive values present
safe_log <- function(x, tiny = 1e-6) {
  x_num <- as.numeric(x)
  if (all(is.na(x_num))) return(x_num)
  if (any(x_num <= 0, na.rm = TRUE)) {
    offset <- abs(min(x_num, na.rm = TRUE)) + tiny
    warning("Non-positive values found; adding offset = ", signif(offset, 6), " before log.")
    x_num <- x_num + offset
  }
  return(log(x_num))
}

# apply log to non-rate variables (excluding date)
for (var in names(df)) {
  if (!(var %in% rate_variables) && var != "date") {
    df[[var]] <- safe_log(df[[var]])
  }
}

# convert to growth rates (percent) for the same non-rate variables
for (var in names(df)) {
  if (!(var %in% rate_variables) && var != "date") {
    df[[var]] <- 100 * (df[[var]] - dplyr::lag(df[[var]], 1))
  }
}

# remove initial rows with NA introduced by growth calculation if needed
# (subsequent code already removes rows with NA in gdp)

#check if stationary using functions from stationary.R
for (var in setdiff(colnames(df), "date")) {  # excluding date column
  cat("Checking stationarity for:", var, "\n")
  ts_data <- ts(df[[var]])
  stationary_ts <- make_stationary(ts_data, use_log = TRUE)   # returns aligned vector (same length)
  df[[var]] <- as.numeric(stationary_ts)                     # assign directly, no extra NA
}
df <- df %>% filter(!is.na(gdp))  # remove first row with NA

# ---------------------- lags ---------------------#

# possible lags: 1, 4 (one year), 8 (two years)
lags <- c(1, 4)


# ------------------ selected variables ------------------#

#selected variables for BVAR model

selected_variables_1 <- c("gdp", "cpi", "wkfreuro", "consp", "consg", 
                          "ifix", "icnstr", "ime", "exc1", "imc1", "ltot", "uroff", "wage", 
                          "srate", "poilusd", "wd", "pcioecd", "vaabcde", "vaghji")

#from this initial set, we will also try a smaller set of variables

selected_variables_2 <- c("gdp", "cpi", "wkfreuro", "consp", "consg", "ifix", 
                          "exc1", "imc1", "ltot", "uroff", "wage", "srate", 
                          "poilusd", "wd", "pcioecd")

selected_variables_3 <- c("gdp", "srate")

#-------------------------------------------------------------------
# List of variable sets to test
variable_sets <- list(
  #set_1 = selected_variables_1,
  #set_2 = selected_variables_2,
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

# Diffuse (flat) prior: approssimazione per OLS usando Minnesota con lambda molto piccolo
diffuse_prior_bv <- bv_priors(
  hyper = c("lambda", "alpha"),
  mn = bv_minnesota(
    lambda = bv_lambda(mode = 1e-6, sd = 1e-6, min = 1e-12, max = 1),
    alpha  = bv_alpha(mode = 4)
  )
)

#priors
mn_prior_bv <- bv_priors(
  hyper = c("lambda", "alpha"),
  mn = mn
)

soc_prior_bv <- bv_priors(
  hyper = c("lambda", "alpha", "psi"),
  mn = mn,
  soc = soc
)


#list of priors
#lista di prior pronta per essere passata a BVAR::bvar
prior_configs <- list(
  mn      = mn_prior_bv,
  soc     = soc_prior_bv,
  diffuse = diffuse_prior_bv
)

rolling_window_size <- 60  # e.g., 60 quarters (15 years)

#---------------------- end of setup ---------------------#
compute_rmse <- function(data, variables, lags, prior, model_type = c("bvar", "var")) {
  model_type <- match.arg(model_type)
  errors_list <- list()
  errors_matrix <- NULL
  valid_windows <- 0
  errors_summary <- character()
  invalid_forecast_vars <- character()
  
  seq_i <- seq(rolling_window_size, nrow(data) - 1)
  for (i in seq_i) {
    train <- data[(i - rolling_window_size + 1):i, variables, drop = FALSE]
    test_row <- data[i + 1, variables, drop = FALSE]
    
    # coerce and validate training & test
    train_num <- as.data.frame(lapply(train, as.numeric))
    if (any(!is.finite(as.matrix(train_num))) || any(is.na(train_num))) {
      errors_summary <- c(errors_summary, sprintf("window %d skipped: NA/non-finite in train", i))
      next
    }
    if (nrow(train_num) <= max(1, lags)) {
      errors_summary <- c(errors_summary, sprintf("window %d skipped: insufficient obs (n=%d)", i, nrow(train_num)))
      next
    }
    test_vals <- as.numeric(test_row)
    if (any(is.na(test_vals)) || any(!is.finite(test_vals))) {
      errors_summary <- c(errors_summary, sprintf("window %d skipped: NA/non-finite in test", i))
      next
    }
    
    # fit & forecast (capture any error)
    if (model_type == "bvar") {
      fit <- tryCatch(BVAR::bvar(data = train_num, lags = lags, priors = prior, ndraw = 1000, burnin = 500),
                      error = function(e) { errors_summary <<- c(errors_summary, paste0("bvar error window ", i, ": ", conditionMessage(e))); NULL })
      if (is.null(fit)) next
      fc <- tryCatch(predict(fit, h = 1), error = function(e) { errors_summary <<- c(errors_summary, paste0("predict(bvar) error window ", i, ": ", conditionMessage(e))); NULL })
      if (is.null(fc) || is.null(fc$fcst)) next
      forecasted_mean <- sapply(fc$fcst, function(el) el$mean[1])
    } else {
      fit <- tryCatch(vars::VAR(train_num, p = max(1, lags), type = "const"),
                      error = function(e) { errors_summary <<- c(errors_summary, paste0("VAR error window ", i, ": ", conditionMessage(e))); NULL })
      if (is.null(fit)) next
      fc <- tryCatch(predict(fit, n.ahead = 1, ci = 0.95),
                     error = function(e) { errors_summary <<- c(errors_summary, paste0("predict(VAR) error window ", i, ": ", conditionMessage(e))); NULL })
      if (is.null(fc) || is.null(fc$fcst)) next
      forecasted_mean <- sapply(variables, function(v) {
        m <- fc$fcst[[v]]
        if (!is.null(m) && "fcst" %in% colnames(m)) return(as.numeric(m[1, "fcst"])) else return(NA_real_)
      })
    }
    
    # validate forecasted values
    if (any(is.na(forecasted_mean)) || any(!is.finite(forecasted_mean))) {
      bad_vars <- variables[which(!is.finite(forecasted_mean) | is.na(forecasted_mean))]
      invalid_forecast_vars <- unique(c(invalid_forecast_vars, bad_vars))
      errors_summary <- c(errors_summary, sprintf("window %d: invalid forecast for vars: %s", i, paste(bad_vars, collapse = ",")))
      next
    }
    
    # accumulate errors
    err <- forecasted_mean - test_vals
    errors_matrix <- if (is.null(errors_matrix)) matrix(err, nrow = 1) else rbind(errors_matrix, err)
    valid_windows <- valid_windows + 1
  }
  
  if (is.null(errors_matrix) || valid_windows == 0) {
    return(list(rmse = NA_real_, valid_windows = valid_windows, invalid_forecast_vars = invalid_forecast_vars, errors_summary = errors_summary))
  }
  rmse_by_var <- sqrt(colMeans(errors_matrix^2))
  overall_rmse <- mean(rmse_by_var)
  return(list(rmse = overall_rmse, rmse_by_var = rmse_by_var, valid_windows = valid_windows, invalid_forecast_vars = invalid_forecast_vars, errors_summary = errors_summary))
}


#–-------------------- BVAR/VAR RMSE computation ---------------------#

results <- list()


for (lag in lags) {
  message("Testing lag: ", lag)
  for (var_set_name in names(variable_sets)) {
    vars <- variable_sets[[var_set_name]]
    for (prior_name in names(prior_configs)) {
      prior_conf <- prior_configs[[prior_name]]
      model_type <- if (prior_name == "diffuse") "var" else "bvar"
      res <- compute_rmse(data = df, variables = vars, lags = lag, prior = prior_conf, model_type = model_type)
      
      # store full result object
      results[[paste(lag, var_set_name, prior_name, sep = "_")]] <- res
      message(sprintf("Lag:%s VarSet:%s Prior:%s Model:%s -> RMSE:%s valid_windows:%d invalid_vars:%s",
                      lag, var_set_name, prior_name, model_type,
                      ifelse(is.na(res$rmse), "NA", format(res$rmse, digits = 6)),
                      res$valid_windows,
                      ifelse(length(res$invalid_forecast_vars)==0, "none", paste(res$invalid_forecast_vars, collapse=","))))
    }
  }
}

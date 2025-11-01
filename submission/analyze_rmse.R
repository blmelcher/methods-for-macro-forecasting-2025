#analyze results
library(dplyr)
library(readr)

# Carica il file RData prodotto da rmse.R (path relativo alla cartella submission)
if (!file.exists("output/results_df.RData")) {
  stop("File 'output/results_df.RData' non trovato. Esegui prima rmse.R per generare 'output/results_df.RData'.")
}
load("output/results_df.RData")  # carica l'oggetto results_df

summary_prior = results_df %>%
  group_by(prior)


summary_set = summary_prior %>%
    group_by(variable_set) %>%
    summarize(mean_rmse = mean(rmse, na.rm = TRUE), n = n(), .groups = "drop")

# calcola la media dell'RMSE per (prior, variable_set)
results_prior_set <- results_df %>%
  group_by(prior, variable_set) %>%
  summarize(mean_rmse = mean(rmse, na.rm = TRUE), n = n(), .groups = "drop")

# stampa tabella riassuntiva (opzionale)
print(results_prior_set)

# find best combination overall
best_prior_set <- results_prior_set %>%
  slice_min(order_by = mean_rmse, n = 1, with_ties = FALSE)

cat("Best combination (prior + variable_set) for mean RMSE:\n")
print(best_prior_set)

# best combination for each prior
best_per_prior <- results_prior_set %>%
  group_by(prior) %>%
  slice_min(order_by = mean_rmse, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Best variable_set for each prior:\n")
print(best_per_prior)



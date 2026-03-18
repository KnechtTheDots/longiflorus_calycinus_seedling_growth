# the first three lines create the log file for the snakemake
# workflow
log <- file(snakemake@log[['log']], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(loo)
library(cmdstanr)


mod1 <- cmdstan_model(snakemake@params[['rgr']])
mod2 <- cmdstan_model(snakemake@params[['rgr_size']])
mod3 <- cmdstan_model(snakemake@params[['rgr_size_seed']])


# read in the data
long <- read.csv(snakemake@input[['long']])
short <- read.csv(snakemake@input[['short']])

print("Data read in!")


# create a list of the data for stan
dat <- list(
  n_long = nrow(long),
  n_short = nrow(short),
  n_pred = 100,
  area = long$seedling_area,
  age = long$age,
  id_long = long$id,
  survive = short$survive,
  id_short = short$id,
  final_size = short$seedling_area,
  seed = short$seed_size,
  age_max = short$age,
  z_tilde = seq(-4,1.5, l = 100)
)

num_chains = ifelse(parallel::detectCores() >= 4, 4, parallel::detectCores())

fit1 <- mod1$sample(data = dat, chains = 4, parallel_chains = num_chains)
fit2 <- mod2$sample(data = dat, chains = 4, parallel_chains = num_chains)
fit3 <- mod3$sample(data = dat, chains = 4, parallel_chains = num_chains)


fit1_loo <- fit1$loo()
fit2_loo <- fit2$loo()
fit3_loo <- fit3$loo()

fit1$diagnostic_summary()
fit2$diagnostic_summary()
fit3$diagnostic_summary()

compared <- loo_compare(fit1_loo, fit2_loo, fit3_loo)
write.csv(compared, snakemake@output[['compare']])

print(compared)
mods <- list(fit1_loo, fit2_loo, fit3_loo)
weights <- loo_model_weights(mods)
write.csv(weights, snakemake@output[['weights']])

print(weights)

on.exit(sink())

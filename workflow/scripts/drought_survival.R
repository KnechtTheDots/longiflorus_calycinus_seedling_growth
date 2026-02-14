# the first three lines create the log file for the snakemake
# workflow
log <- file(snakemake@log[['log']], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# read in the data
long <- read.csv(snakemake@input[['long']])
short <- read.csv(snakemake@input[['short']])

long <- read.csv("results/wrangled/drought_long.csv")
short <- read.csv("results/wrangled/drought_short.csv")

print("Data read in!")

# create a list of the data for stan
dat <- list(
  n_long = nrow(long),
  n_short = nrow(short),
  area = long$seedling_area,
  age = long$age,
  id_long = long$id,
  survive = short$survive,
  id_short = short$id,
  final_size = short$seedling_area,
  seed = short$seed_size,
  age_max = short$age
)

# compile the stan file into the executable 
mod <- cmdstanr::cmdstan_model(snakemake@params[['stan']])

num_cores <- ifelse(parallel::detectCores() >= 4, 4, parallel::detectCores())


# sample from the posterior
fit <- mod$sample(
  data = dat,
  chains = 4,
  parallel_chains = num_cores
)


# save the fit object for later analysis
fit$save_object(snakemake@output[['fit']])

print("Fit object written to file")

on.exit(sink())
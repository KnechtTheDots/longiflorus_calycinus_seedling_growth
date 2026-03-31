# the first three lines create the log file for the snakemake
# workflow
log <- file(snakemake@log[['log']], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(cmdstanr)

mod <- cmdstan_model(snakemake@input[['sze']])
mod <- cmdstan_model(snakemake@input[['rgr']])
mod <- cmdstan_model(snakemake@input[['rgr_size']])

write.csv("models compiled!", snakemake@output[[1]])





on.exit(sink())
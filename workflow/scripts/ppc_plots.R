# the first three lines create the log file for the snakemake
# workflow
log <- file(snakemake@log[['log']], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# load in libraries
library(ggplot2)
library(dplyr)
library(bayesplot)
library(posterior)
library(cmdstanr)
library(cowplot)

print("Packages read in")

# read in the cmdstanfit object
fit <- readRDS(snakemake@input[['fit']])
# read in the short dataset
short <- read.csv(snakemake@input[['short']])

print("Data read in")

# get the posterior predictions of the final sizes
y_pred <- fit$draws("y_rep", format = "df")[,1:nrow(short)]

print("starting calcs")

# calculate the mean and, as well as the upper and lower 90% credible interval
# of the posterior predictive distribution
mu <- as.numeric(apply(y_pred, 2, mean))
upr <- as.numeric(apply(y_pred, 2, function(x) quantile(x, .95)))
lwr <- as.numeric(apply(y_pred, 2, function(x) quantile(x, .05)))

# plot the observed log sizes vs the predicted log sizes
obs_v_pred <- data.frame(mu, upr, lwr, obs = log(short$seedling_area)) %>% 
  ggplot(aes(x = obs, y = mu)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_errorbar(aes(x = obs, ymax = upr, ymin = lwr)) +
  labs(x = "Observed size",
       y = "Predicted size") +
  theme_classic()


print("Data prediction plot made")

y <- log(short$seedling_area)
# get the posterior predictive, but use matrix this time for the bayesplot functions
y_rep <- fit$draws("y_rep", format = "matrix")

# sample 100 random sets of draws for the density overlay plot
samps <- sample(1:nrow(y_rep), 100)

# plot density of the predictions vs the actual
dens_plot <- ppc_dens_overlay(y, yrep = y_rep[samps,])
# plot histogram of standard deviation of posterior predictions vs observed
sd_hist <- ppc_stat(y, yrep = y_rep, stat = "sd")
# plot histogram of mean of posterior predictions vs observed
mu_hist <- ppc_stat(y, yrep = y_rep, stat = "mean")
# plot histogram of minimum of posterior predictions vs observed
min_hist <- ppc_stat(y, yrep = y_rep, stat = "min")
# plot histogram of maximum of posterior predictions vs observed
max_hist <- ppc_stat(y, yrep = y_rep, stat = "max")

# read in the summary of fit to gett effective sample sizes and rhat
d <- fit$summary()

# plot bulk ess vs rhat
rhat_bulk <- d %>% 
  ggplot(aes(x = ess_bulk, y = rhat)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Bulk ESS",
       y = "Rhat") +
  geom_hline(yintercept = 1.01, color = "red") +
  geom_vline(xintercept = 400, color = "red")



# plot tail ess bs rhat
rhat_tail <- d %>% 
  ggplot(aes(x = ess_tail, y = rhat)) + 
  geom_point() +
  theme_minimal() +
  labs(x = "Tail ESS",
       y = "Rhat") +
  geom_hline(yintercept = 1.01, color = "red") +
  geom_vline(xintercept = 400, color = "red")

grid <- plot_grid(obs_v_pred, dens_plot, mu_hist, sd_hist,
          min_hist, max_hist, rhat_bulk, rhat_tail,
          ncol = 4, labels = "AUTO")

ggsave(snakemake@output[[1]], grid, device = "svg", width = 12, height = 8)


on.exit(sink())
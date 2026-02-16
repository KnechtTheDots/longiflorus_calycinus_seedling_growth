# the first three lines create the log file for the snakemake
# workflow
log <- file(snakemake@log[['log']], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(ggplot2)
library(dplyr)
library(cmdstanr)

this_theme <- theme(axis.title = element_text(size = 15, face = "bold"),
                    axis.text = element_text(face = "bold"))

get_quant <- function(x, quant){
  as.numeric(apply(x, 2, function(z) quantile(z, quant)))
}

# read in the fitted stan model and short dataset
fit <- readRDS(snakemake@input[['fit']])
short <- read.csv(snakemake@input[['short']])

print("Data read in.")

# standardize the final size
short <- short %>% 
  mutate(z_size = (seedling_area - mean(seedling_area))/sd(seedling_area))

# create the predictor variables that the stan model predicted over
size_pred = seq(-2,2.5,l=100)
rgr_pred = seq(.125, .325, l = 100)

# extract the posterior draws of survival probability across
p_size_pred <- fit$draws("p_size_pred", format = "df")[,1:length(size_pred)]

# get the mean survival probability at each size
mu <- as.numeric(apply(p_size_pred, 2, mean))
# get the 97.5% prob at each size
upr <- get_quant(p_size_pred, .975)
# get the 2.5% prob at each size
lwr <- get_quant(p_size_pred, .025)
# get the 75% prob at each size
upr.5 <- get_quant(p_size_pred, .75)
# get the 25% prob at each size
lwr.5 <- get_quant(p_size_pred, .25)

# plot survival vs size
data.frame(mu, upr, lwr, upr.5, lwr.5, size = size_pred) %>% 
  ggplot(aes(x = size, y = mu)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(x = size, ymin = lwr, ymax = upr), alpha = .25) +
  geom_ribbon(aes(x = size, ymin = lwr.5, ymax = upr.5), alpha = .25) +
  geom_jitter(data = short, aes(x = z_size, y = survive), width = 0, height = .01) +
  labs(x = "Standardized Size",
       y = "Probability of Survival") +
  theme_minimal() +
  this_theme

# save the plot
ggsave(snakemake@output[[1]], device = "svg", width = 12, height = 8)

# extract the relative growth rate of each individual and calculate the mean
rgr <- fit$draws("rgr", format = "df")[,1:nrow(short)]
mu_rgr <- as.numeric(apply(rgr, 2, mean))

# extract the predicted survival across relative growth rates at day 14
p_rgr <- fit$draws("p_rgr_age14_pred", format = "df")[,1:length(rgr_pred)]

# create a data frame to include the points for each individuals relative growth rate 
# on the plot and color them by age
rgr_surv <- data.frame(rgr = mu_rgr, survive = short$survive, age = factor(short$age, levels = c(11,12,13,14)))

# get the mean and 95%, 50% intervals
mu <- as.numeric(apply(p_rgr, 2, mean))
upr <- get_quant(p_rgr, .975)
lwr <- get_quant(p_rgr, .025)
upr.5 <- get_quant(p_rgr, .75)
lwr.5 <- get_quant(p_rgr, .25)


# plot the survival vs relative growth
data.frame(mu, upr, lwr, upr.5, lwr.5, rgr_pred) %>% 
  ggplot(aes(x = rgr_pred, y = mu)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(x = rgr_pred, ymax = upr, ymin = lwr), alpha = .25) +
  geom_ribbon(aes(x = rgr_pred, ymax = upr.5, ymin = lwr.5), alpha = .25) +
  geom_jitter(data = rgr_surv, aes(x = rgr, y = survive, color = age), width = 0, height = .02) +
  scale_colour_viridis_d() +
  labs(x = "Relative Growth Rate",
       y = "Probability of Survival",
       color = "Age") +
  theme_minimal() +
  this_theme

ggsave(snakemake@output[[2]], device = "svg", width = 12, height = 8)

# compare suvival vs growth rate at different ages
p_rgr_11 <- fit$draws("p_rgr_age11_pred", format = "df")[,1:length(rgr_pred)]
p_rgr_12 <- fit$draws("p_rgr_age12_pred", format = "df")[,1:length(rgr_pred)]
p_rgr_13 <- fit$draws("p_rgr_age13_pred", format = "df")[,1:length(rgr_pred)]
p_rgr_14 <- fit$draws("p_rgr_age14_pred", format = "df")[,1:length(rgr_pred)]
mu_11 <- as.numeric(apply(p_rgr_11, 2, mean))
mu_12 <- as.numeric(apply(p_rgr_12, 2, mean))
mu_13 <- as.numeric(apply(p_rgr_13, 2, mean))
mu_14 <- as.numeric(apply(p_rgr_14, 2, mean))

# plot survival vs age vs rgr
data.frame(mu_11, mu_12, mu_13, mu_14, rgr_pred) %>% 
  tidyr::pivot_longer(1:4, names_to = "age", values_to = "p") %>% 
  mutate(age = factor(gsub("mu_", "", age), levels = c(11,12,13,14))) %>% 
  ggplot(aes(x = rgr_pred, y = p, color = age)) +
  geom_line(linewidth = 1) +
  scale_colour_viridis_d() +
  labs(x = "Relative Growth Rate",
       y = "Probability of Survival",
       color = "Age") +
  theme_minimal() +
  this_theme

ggsave(snakemake@output[[3]], device = "svg", width = 12, height = 8)

print("All plots made")

on.exit(sink())




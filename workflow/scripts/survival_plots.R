library(ggplot2)
library(dplyr)
library(cmdstanr)

get_quant <- function(x, quant){
  as.numeric(apply(x, 2, function(z) quantile(z, quant)))
}

fit <- readRDS("results/drought_survival.RDS")

fit$diagnostic_summary()

short <- read.csv("results/wrangled/drought_short.csv")

short <- short %>% 
  mutate(z_size = (seedling_area - mean(seedling_area))/sd(seedling_area))

size_pred = seq(-2,2.5,l=100)
rgr_pred = seq(.125, .325, l = 100)

p_size_pred <- fit$draws("p_size_pred", format = "df")[,1:length(size_pred)]

mu <- as.numeric(apply(p_size_pred, 2, mean))
upr <- get_quant(p_size_pred, .975)
lwr <- get_quant(p_size_pred, .025)
upr.5 <- get_quant(p_size_pred, .75)
lwr.5 <- get_quant(p_size_pred, .25)

data.frame(mu, upr, lwr, upr.5, lwr.5, size = size_pred) %>% 
  ggplot(aes(x = size, y = mu)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(x = size, ymin = lwr, ymax = upr), alpha = .25) +
  geom_ribbon(aes(x = size, ymin = lwr.5, ymax = upr.5), alpha = .25) +
  geom_point(data = short, aes(x = z_size, y = survive)) +
  theme_minimal()


p_rgr <- fit$draws("p_rgr_age14_pred", format = "df")[,1:length(rgr_pred)]

rgr <- fit$draws("rgr", format = "df")[,1:nrow(short)]

mu_rgr <- as.numeric(apply(rgr, 2, mean))

rgr_surv <- data.frame(rgr = mu_rgr, survive = short$survive, age = factor(short$age, levels = c(11,12,13,14)))

mu <- as.numeric(apply(p_rgr, 2, mean))
upr <- get_quant(p_rgr, .975)
lwr <- get_quant(p_rgr, .025)
upr.5 <- get_quant(p_rgr, .75)
lwr.5 <- get_quant(p_rgr, .25)



data.frame(mu, upr, lwr, upr.5, lwr.5, rgr_pred, mu_10, mu_12) %>% 
  ggplot(aes(x = rgr_pred, y = mu)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(x = rgr_pred, ymax = upr, ymin = lwr), alpha = .25) +
  geom_ribbon(aes(x = rgr_pred, ymax = upr.5, ymin = lwr.5), alpha = .25) +
  geom_jitter(data = rgr_surv, aes(x = rgr, y = survive, color = age), width = 0, height = .02) +
  scale_colour_viridis_d() +
  theme_minimal()

p_rgr_11 <- fit$draws("p_rgr_age11_pred", format = "df")[,1:length(rgr_pred)]
p_rgr_12 <- fit$draws("p_rgr_age12_pred", format = "df")[,1:length(rgr_pred)]
p_rgr_13 <- fit$draws("p_rgr_age13_pred", format = "df")[,1:length(rgr_pred)]
p_rgr_14 <- fit$draws("p_rgr_age14_pred", format = "df")[,1:length(rgr_pred)]
mu_11 <- as.numeric(apply(p_rgr_11, 2, mean))
mu_12 <- as.numeric(apply(p_rgr_12, 2, mean))
mu_13 <- as.numeric(apply(p_rgr_13, 2, mean))
mu_14 <- as.numeric(apply(p_rgr_14, 2, mean))

data.frame(mu_11, mu_12, mu_13, mu_14, rgr_pred) %>% 
  tidyr::pivot_longer(1:4, names_to = "age", values_to = "p") %>% 
  mutate(age = factor(gsub("mu_", "", age), levels = c(11,12,13,14))) %>% 
  ggplot(aes(x = rgr_pred, y = p, color = age)) +
  geom_line(linewidth = 1) +
  scale_colour_viridis_d() +
  theme_minimal()








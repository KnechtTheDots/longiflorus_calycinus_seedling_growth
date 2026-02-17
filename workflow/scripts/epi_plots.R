# the first three lines create the log file for the snakemake
# workflow
log <- file(snakemake@log[['log']], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggbeeswarm)
library(cowplot)

this_theme <- theme(axis.title = element_text(size = 10, face = "bold"),
                    axis.text = element_text(face = "bold"))



plot_epi_test <- function(df, trait_id){
  df %>% 
    pivot_longer(1:5, names_to = "line", values_to = "trait") %>% 
    group_by(line) %>% 
    summarise(mu = mean(trait),
              upr = quantile(trait, .975),
              lwr = quantile(trait, .025),
              upr.5 = quantile(trait, .75),
              lwr.5 = quantile(trait, .25)) %>% 
    mutate(group = c(1,2,3,4,1),
           prop_lon = c(0, .5, .48, .52, 1)) %>% 
    ggplot(aes(x = prop_lon, y = mu, color = line, group = group)) +
    geom_line(color = "black", linewidth = 1) +
    geom_errorbar(aes(x = prop_lon, ymax = upr, ymin = lwr), linewidth = 1, width = 0) +
    geom_errorbar(aes(x = prop_lon, ymax = upr.5, ymin = lwr.5),
                  linewidth = 1.5, width= 0) +
    labs(x = "Proportion Longiflorus",
         y = trait_id) +
    scale_color_manual(values = c("black", "red", "black", "blue", "black")) +
    theme_minimal() +
    theme(legend.position = "none") +
    this_theme
}


rgr <- read.csv(snakemake@input[['rgr']])

p_rgr <- plot_epi_test(rgr, trait_id = "Relative Growth Rate")

day_4 <- read.csv(snakemake@input[['day_4']])

p_d4 <- plot_epi_test(day_4, trait_id = "Area Day 4")

day_17 <- read.csv(snakemake@input[['day_17']])

p_d17 <- plot_epi_test(day_17, trait_id = "Area Day 17")

surv <- read.csv(snakemake@input[['survive']])

p_surv <- plot_epi_test(surv, trait_id = "Proportion Survive")

height <- read.csv(snakemake@input[['height']])

p_height <- plot_epi_test(height, trait_id = "Height")


delta <- function(f2, f1, p1, p2){
  return(f2 - (.5*f1 + .25*(p1 + p2)))
}



delta_f2 <- data.frame(Survive = delta(surv$f2, surv$f1, surv$lon, surv$cal)/surv$e_f2,
           day_4 = delta(day_4$f2, day_4$f1, day_4$lon, day_4$cal)/day_4$e_f2,
           day_17 = delta(day_17$f2, day_17$f1, day_17$lon, day_17$cal)/day_17$e_f2,
           rgr = delta(rgr$f2, rgr$f1, rgr$lon, rgr$cal)/rgr$e_f2,
           height = delta(height$f2, height$f1, height$lon, height$cal)/height$e_f2) %>% 
  pivot_longer(1:5, names_to = "trait", values_to = "delta") %>% 
  group_by(trait) %>% 
  summarise(mu = mean(delta),
            upr = quantile(delta, .975),
            lwr = quantile(delta, .025),
            upr.5 = quantile(delta, .75),
            lwr.5 = quantile(delta, .25)) %>% 
  mutate(Trait = c("Survive", "Day 17", "Day 4", "Height", "RGR"),
         Trait = factor(Trait, levels = c("Survive", "Day 4", "Day 17", "RGR", "Height"))) %>% 
  ggplot(aes(x = Trait, y = mu)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(x = Trait, ymax = upr, ymin = lwr), width = 0, linewidth = 1) +
  geom_errorbar(aes(x = Trait, ymax = upr.5, ymin = lwr.5), width = 0,
                linewidth = 1.5) +
  labs(y = "Deviation/E(F2)",
       x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(face = "bold"))



plot_grid(p_surv, p_d4, p_d17, p_rgr, p_height, delta_f2, ncol = 3, labels = "AUTO")

ggsave(snakemake@output[['epi']], device = "svg", width = 12, height = 8)


d <- read.csv(snakemake@input[['traits']])

d$rgr <- (log(d$day_17) - log(d$day_4))/(17-4)

d4 <- d %>% 
  ggplot(aes(x = line, y = day_4)) +
  geom_beeswarm() +
  theme_minimal() +
  labs(x = "",
       y = "Area Day 4") +
  this_theme +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

d17 <- d %>% 
  ggplot(aes(x = line, y = day_17)) +
  geom_beeswarm() +
  theme_minimal() +
  labs(x = "",
       y = "Area Day 17") +
  this_theme +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

drgr <- d %>% 
  ggplot(aes(x = line, y = rgr)) +
  geom_beeswarm() +
  theme_minimal() +
  labs(x = "",
       y = "Relative Growth Rate") +
  this_theme +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

dheight <- d %>% 
  ggplot(aes(x = line, y = height)) +
  geom_beeswarm() +
  theme_minimal() +
  labs(x = "",
       y = "Height") +
  this_theme +
  theme(axis.text.x = element_text(size = 12, face = "bold"))


plot_grid(d4, d17, drgr, dheight, ncol = 2, labels = "AUTO")
ggsave(snakemake@output[['obs']], device = "svg", width = 12, height = 12)

on.exit(sink())
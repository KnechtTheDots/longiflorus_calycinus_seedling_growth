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

this_theme <- theme(axis.title = element_text(size = 13, face = "bold"),
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
           prop_lon = c(0, .5, .47, .53, 1)) %>% 
    ggplot(aes(x = prop_lon, y = mu, color = line, group = group)) +
    geom_line(color = "black", linewidth = 1) +
    geom_errorbar(aes(x = prop_lon, ymax = upr, ymin = lwr), linewidth = 1, width = 0) +
    geom_errorbar(aes(x = prop_lon, ymax = upr.5, ymin = lwr.5),
                  linewidth = 1.75, width= 0) +
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
            lwr.5 = quantile(delta, .25),
            prop = round(mean(delta < 0), 2)) %>%
  ungroup() %>% 
  mutate(prop_place = max(upr)*1.2) %>% 
  mutate(Trait = c("Survive", "Day 17", "Day 4", "Height", "RGR"),
         Trait = factor(Trait, levels = c("Survive", "Day 4", "Day 17", "RGR", "Height"))) %>% 
  ggplot(aes(x = Trait, y = mu)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_errorbar(aes(x = Trait, ymax = upr, ymin = lwr), width = 0, linewidth = 1) +
  geom_errorbar(aes(x = Trait, ymax = upr.5, ymin = lwr.5), width = 0,
                linewidth = 1.75) +
  geom_text(aes(x = Trait, y = prop_place, label = prop),
            fontface = "bold") +
  labs(y = expression(Delta/E(F[2])),
       x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold", vjust = -3),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"))



grid <- plot_grid(p_surv, p_d4, p_d17, p_rgr, p_height, delta_f2, ncol = 3, labels = "AUTO")

ggsave(snakemake@output[['epi']], grid, device = "jpg", width = 12, height = 8)


this_theme <- theme(axis.title = element_text(size = 13, face = "bold"),
                    axis.text = element_text(face = "bold"))

d <- read.csv(snakemake@input[['traits']])

d$rgr <- (log(d$day_17) - log(d$day_4))/(17-4)

# get sample sizes for day for to plot
n4 <- d %>% 
  drop_na(day_4) %>% 
  group_by(line) %>% 
  summarise(n = n()) %>% 
  ungroup()

# height of the sample size labels
n4$lab_y <- max(d$day_4, na.rm = T) + .25

d4 <- d %>% 
  ggplot(aes(x = line, y = day_4)) +
  geom_beeswarm() +
  theme_minimal() +
  geom_text(data = n4, aes(x = line, y = lab_y, label = n),
            fontface = "bold") +
  labs(x = "",
       y = "Area Day 4") +
  this_theme +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

# get sample sizes for day for to plot
n17 <- d %>%
  drop_na(day_17) %>%
  group_by(line) %>%
  summarise(n = n()) %>%
  ungroup()

# height of the sample size labels
n17$lab_y <- max(d$day_17, na.rm = T) + 5

d17 <- d %>% 
  ggplot(aes(x = line, y = day_17)) +
  geom_beeswarm() +
  theme_minimal() +
  geom_text(data = n17, aes(x = line, y = lab_y, label = n),
            fontface = "bold") +
  labs(x = "",
       y = "Area Day 17") +
  this_theme +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

# get sample sizes to plot
nrgr <- d %>%
    drop_na(rgr) %>%
    group_by(line) %>%
    summarise(n = n()) %>%
    ungroup()

# height of the sample size labels
nrgr$lab_y <- max(d$rgr, na.rm = T) + .025
  
drgr <- d %>% 
  ggplot(aes(x = line, y = rgr)) +
  geom_beeswarm() +
  theme_minimal() +
  geom_text(data = nrgr, aes(x = line, y = lab_y, label = n),
            fontface = "bold") +
  labs(x = "",
       y = "Relative Growth Rate") +
  this_theme +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

# get sample sizes for day for to plot
nheight <- d %>%
  drop_na(height) %>%
  group_by(line) %>%
  summarise(n = n()) %>%
  ungroup()

# height of the sample size labels
nheight$lab_y <- max(d$height, na.rm = T) + 2.5

dheight <- d %>% 
  ggplot(aes(x = line, y = height)) +
  geom_beeswarm() +
  theme_minimal() +
  geom_text(data = nheight, aes(x = line, y = lab_y, label = n),
            fontface = "bold") +
  labs(x = "",
       y = "Height") +
  this_theme +
  theme(axis.text.x = element_text(size = 12, face = "bold"))


grid <- plot_grid(d4, d17, drgr, dheight, ncol = 2, labels = "AUTO")
ggsave(snakemake@output[['obs']], grid, device = "jpg", width = 12, height = 8)


this_theme <- theme(axis.title = element_text(size = 12, face = "bold"),
                    axis.text = element_text(face = "bold"))

mean_diffs <- function(df, trait){
  data.frame("F2-F1" = df$f2 - df$f1,
             "F2-LON" = df$f2 - df$lon,
             "F2-CAL" = df$f2 - df$cal,
             "F1-LON" = df$f1 - df$lon,
             "F1-CAL" = df$f1 - df$cal,
             "LON-CAL" = df$lon - df$cal) %>% 
    pivot_longer(1:6, names_to = "comp", values_to = "diff") %>% 
    group_by(comp) %>% 
    mutate(comp = case_when(comp == "F2.F1" ~ "F2-F1",
                            comp == "F2.LON" ~ "F2-LON",
                            comp == "F2.CAL" ~ "F2-CAL",
                            comp == "F1.LON" ~ "F1-LON",
                            comp == "F1.CAL" ~ "F1-CAL",
                            comp == "LON.CAL" ~ "LON-CAL"),
           comp = factor(comp, levels = c("F2-F1","F2-LON","F2-CAL",
                                          "F1-LON","F1-CAL","LON-CAL"))) %>% 
    summarise(mu = mean(diff),
              upr = quantile(diff, .975),
              lwr = quantile(diff, .025),
              upr.5 = quantile(diff, .75),
              lwr.5 = quantile(diff, .25),
              prop = round(mean(diff < 0), 2),
              prop_place = max(upr)*1.2) %>%
    ungroup() %>% 
    mutate(prop_place = max(prop_place)) %>% 
    ggplot(y = comp, x = mu) +
    geom_vline(xintercept = 0, color = "grey") +
    geom_errorbarh(aes(y = comp, xmax = upr, xmin = lwr), height = 0, linewidth = 1) +
    geom_errorbarh(aes(y = comp, xmax = upr.5, xmin = lwr.5), height = 0, linewidth = 1.75) +
    geom_text(aes(x = prop_place, y = comp, label = prop),
              fontface = "bold") +
    labs(y = "",
         x = paste0("Difference in ", trait, " Means")) +
    theme_minimal() +
    this_theme
}


diff_surv <- mean_diffs(surv, "Survival")
diff_d4 <- mean_diffs(day_4, "Size on Day 4")
diff_d17 <- mean_diffs(day_17, "Size on Day 17")
diff_rgr <- mean_diffs(rgr, "RGR")
diff_height <- mean_diffs(height, "Height on Day 122")



#grid <- plot_grid(diff_surv, diff_d4, diff_d17, diff_rgr, diff_height, labels = "AUTO", ncol = 3)

row1 <- plot_grid(diff_surv, diff_d4, diff_d17, ncol = 3, labels = c("A", "B", "C"))

blank <- ggdraw() + theme_void()

row2 <- plot_grid(blank, diff_rgr, diff_height, blank, ncol = 4, rel_widths = c(1,2,2,1), labels = c("", "D", "E", ""))

grid <- plot_grid(row1, row2, ncol = 1)

# with the jpg file type the background renders to black, so need to make it white

final <- ggdraw(grid) +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(snakemake@output[['mean_diffs']], final, device = "jpg", width = 12, height = 8)



on.exit(sink())
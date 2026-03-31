# the first three lines create the log file for the snakemake
# workflow
log <- file(snakemake@log[['log']], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(ggplot2)
library(dplyr)
library(tidyr)
library(cmdstanr)
library(cowplot)

# function for calculating the path coefficients
get_cors <- function(r_rgr_size, r_age_size, r_seed_size, d, e, f){
  c <- ((r_seed_size - f*r_rgr_size) - (e - d*f)*(r_age_size - d*r_rgr_size)/(1-d^2))/
    (1 - f^2 - (e-d*f)^2/(1-d^2))
  b <- (r_age_size - d*r_rgr_size)/(1-d^2) - c*(e-d*f)/(1-d^2)
  a <- r_rgr_size - b*d - c*f
  return(data.frame(a, b, c, d, e, f))
}

this_theme <- theme(axis.title = element_text(size = 12, face = "bold"),
                    axis.text = element_text(face = "bold"))

get_quant <- function(x, quant){
  as.numeric(apply(x, 2, function(z) quantile(z, quant)))
}

# read in the fitted stan model and short dataset
short <- read.csv(snakemake@input[['short']])
long <- read.csv(snakemake@input[['long']])


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

plots <- path_coefs <- list()
for(i in 1:3){
  fit <- readRDS(snakemake@input[[i]])
  
  
  print("Data read in.")
  
  z_tilde <- seq(-3.5,1,l=100)
  
  p <- fit$draws("p_pred", format = "df")[,1:length(z_tilde)]
  
  colnames(p) <- z_tilde
  
  rgr_est <- fit$draws("z_rgr", format = "df")[,1:nrow(short)]
  
  mu_rgr <- apply(rgr_est, 2, mean)
  
  rgr_surv <- data.frame(z_rgr = mu_rgr, survive = short$survive)
  
  plots[[i]] <- p %>% 
    pivot_longer(1:ncol(p), names_to = "z_rgr", values_to = "p") %>% 
    group_by(z_rgr) %>% 
    summarise(mu = mean(p),
              upr = quantile(p, .975),
              lwr = quantile(p, .025),
              upr.5 = quantile(p, .75),
              lwr.5 = quantile(p, .25)) %>% 
    mutate(z_rgr = as.numeric(z_rgr)) %>% 
    ggplot(aes(x = z_rgr, y = mu)) +
    geom_line(linewidth = .75) +
    geom_ribbon(aes(x = z_rgr, ymax = upr, ymin = lwr), alpha = .25) +
    geom_ribbon(aes(x = z_rgr, ymax = upr.5, ymin = lwr.5), alpha = .25) +
    geom_jitter(data = rgr_surv, aes(x = z_rgr, y = survive), width = 0, height = .01) +
    labs(x = "Standardized RGR",
         y = "Survival Probability") +
    theme_minimal() +
    this_theme
  
  
  ## plot the correlations
  
  
  d <- fit$draws("R_phenos", format = "df")$`R_phenos[2,1]`
  e <- fit$draws("R_phenos", format = "df")$`R_phenos[2,3]`
  f <- fit$draws("R_phenos", format = "df")$`R_phenos[1,3]`
  
  rgr <- fit$draws("rgr", format = "df")[,1:nrow(short)]
  size <- fit$draws("final_rep", format = "df")[,1:nrow(short)]
  age <- dat$age_max
  seed <- dat$seed
  
  w <- fit$draws("w", format = "df")[,1:nrow(short)]
  sigma_w <- fit$draws("sigma_w", format = "df")$sigma_w
  
  r_rgr_size <- r_age_size <- r_seed_size <- r_rgr_w <- r_size_w <- c()
  for(j in 1:nrow(rgr)){
    r_rgr_size[j] <- cor(as.numeric(rgr[j,]), as.numeric(size[j,]))
    r_age_size[j] <- cor(age, as.numeric(size[j,]))
    r_seed_size[j] <- cor(seed, as.numeric(size[j,]))
    r_rgr_w[j] <- cor(as.numeric(rgr[j,]), as.integer(w[j,]))
    r_size_w[j] <- cor(as.numeric(size[j,]), as.integer(w[j,]))
  }
  
  cors <- get_cors(r_rgr_size, r_age_size, r_seed_size, d, e, f)
  
  if(i==1){
    cors$g <- r_size_w
    cors$h <- 0
  }
  
  if(i==2){
    cors$g <- 0
    cors$h <- r_rgr_w
  }
  
  if(i==3){
    cors$g <- (r_size_w - r_rgr_size*r_rgr_w)/(1 - r_rgr_size^2)
    
    cors$h <- r_rgr_w - cors$g*r_rgr_size
  }
  
  # multiply the isolated correlations by the standard deviation for fitness
  # to get the standardized selection differentials
  selection_diffs <- data.frame(
    RGR = (cors$h + cors$a * cors$g) * sigma_w,
    Age = (cors$b * cors$g) * sigma_w,
    SeedSize = (cors$c * cors$g) * sigma_w
  )
  
  path_coefs[[i]] <- selection_diffs %>% 
    pivot_longer(1:ncol(selection_diffs), names_to = "path",
                 values_to = "S") %>% 
    group_by(path) %>% 
    summarise(upr = quantile(S, .975),
              lwr = quantile(S, .025),
              upr.5 = quantile(S, .75),
              lwr.5 = quantile(S, .25)) %>%  
    mutate(path = factor(path, levels = c("RGR", "Age", "SeedSize"))) %>% 
    ggplot(aes()) +
    geom_hline(yintercept = 0, color = "grey") +
    geom_errorbar(aes(x = path, ymin = lwr, ymax = upr),
                   width = 0, linewidth = .75) +
    geom_errorbar(aes(x = path, ymin = lwr.5, ymax = upr.5),
                   width = 0, linewidth = 1.5) +
    labs(y = "Selection Differential") +
    ylim(-.12, .6) +
    theme_minimal() +
    theme(axis.title.x = element_blank()) +
    this_theme
  
  

  
}

print("Main loop finished")
  
mods <- list()
for(i in 1:3){
  to_survive <- data.frame(x = c(-1, 0), xend = c(-.06,0),
                           y = c(-.9,.1), yend = c(.9,.9))
  to_size <- data.frame(x = c(-1,0,1), xend = c(-.1,0,.1),
                        y = c(-.9,-.9,-.9), yend = c(-.1,-.1,-.1))
  to_phenos <- data.frame(x = c(-.9,-.9,.1), xend = c(-.1,.9,.9),
                          y = c(-1.1,-1.1,-1.1), yend = c(-1.1,-1.1,-1.1))
  
  paths <- data.frame(variables = letters[1:8],
                      x = c(-.5, .05, .5, -.5, .5, 0, .05, -.5),
                      y = c(-.55, -.55, -.55, -1.1, -1.1, -1.4, .3, .3))
  
  if(i==1){
    to_survive <- to_survive[2,]
    paths <- paths[paths$variables != "h",]
  }
  
  if(i==2){
    to_survive <- to_survive[1,]
    paths <- paths[paths$variables != "g",]
  }
  
  mods[[i]] <- data.frame(variable = c("w", "Size", "RGR", "Age","SeedSize"),
             x = c(0, 0, -1, 0, 1),
             y = c(1, 0, -1, -1, -1)) %>% 
    ggplot(aes(x = x, y = y, label = variable)) +
    geom_text(fontface = "bold") +
    geom_segment(data = to_survive, mapping = aes(x = x, y = y,
                                                  xend = xend, yend = yend),
                 linewidth = .75,
                 inherit.aes = F,
                 arrow = arrow(length = unit(0.03, "npc"))) +
    geom_segment(data = to_size, mapping = aes(x = x,y = y,
                                               xend = xend, yend = yend),
                 linewidth = .75,
                 inherit.aes = F,
                 arrow = arrow(length = unit(0.03, "npc"))) +
    geom_curve(data = to_phenos, mapping = aes(x = x, y = y,
                                               xend = xend, yend = yend),
               curvature = .3,
               linewidth = .75,
               arrow = arrow(length = unit(0.03, "npc"), ends = "both"),
               inherit.aes = F) +
    geom_text(data = paths, mapping = aes(x = x, y = y, label = variables), inherit.aes = F) +
    ylim(-1.5,1) + 
    xlim(-1.05,1.1) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank())

}

print("Path plots made")


# plot size, rgr, size-rgr
final_plot <- plot_grid(plots[[1]], plots[[2]], plots[[3]], 
          mods[[1]], mods[[2]], mods[[3]],
          path_coefs[[1]], path_coefs[[2]], path_coefs[[3]],ncol = 3, labels = "AUTO")

print("final plot made")

ggsave(snakemake@output[[1]], final_plot, device = "svg", width = 12, height = 8)


print("All plots made")

on.exit(sink())




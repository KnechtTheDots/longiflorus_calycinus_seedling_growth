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

this_theme <- theme(axis.title = element_text(size = 15, face = "bold"),
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
    geom_line(linewidth = 1) +
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
  
  survive <- fit$draws("surv_rep", format = "df")[,1:nrow(short)]
  
  r_rgr_size <- r_age_size <- r_seed_size <- r_rgr_surv <- r_size_surv <- c()
  for(j in 1:nrow(rgr)){
    r_rgr_size[j] <- cor(as.numeric(rgr[j,]), as.numeric(size[j,]))
    r_age_size[j] <- cor(age, as.numeric(size[j,]))
    r_seed_size[j] <- cor(seed, as.numeric(size[j,]))
    r_rgr_surv[j] <- cor(as.numeric(rgr[j,]), as.integer(survive[j,]))
    r_size_surv[j] <- cor(as.numeric(size[j,]), as.integer(survive[j,]))
  }
  
  cors <- get_cors(r_rgr_size, r_age_size, r_seed_size, d, e, f)
  
  if(i==1){
    cors$g <- r_size_surv
  }
  
  if(i==2){
    cors$g <- (r_size_surv - r_rgr_size*r_rgr_surv)/(1 - r_rgr_size^2)
    
    cors$h <- r_rgr_surv - cors$g*r_rgr_size
  }
  
  levs <- rev(letters[1:ncol(cors)])
  
  path_coefs[[i]] <- cors %>% 
    pivot_longer(1:ncol(cors), names_to = "path", values_to = "cor") %>% 
    group_by(path) %>% 
    summarise(upr = quantile(cor, .975),
              lwr = quantile(cor, .025),
              upr.5 = quantile(cor, .75),
              lwr.5 = quantile(cor, .25)) %>% 
    mutate(path = factor(path, levels = levs)) %>% 
    ggplot(aes()) +
    geom_vline(xintercept = 0, color = "grey") +
    geom_errorbarh(aes(y = path, xmin = lwr, xmax = upr),
                   height = 0, linewidth = .75) +
    geom_errorbarh(aes(y = path, xmin = lwr.5, xmax = upr.5),
                   height = 0, linewidth = 1.5) +
    theme_minimal() +
    theme(axis.title.y = element_blank())
  
}


  
mods <- list()
for(i in 1:3){
  to_survive <- data.frame(x = c(-1, 0), xend = c(-.06,0),
                           y = c(-.94,.06), yend = c(.95,.95))
  to_size <- data.frame(x = c(-1,0,1), xend = c(-.06,0,.06),
                        y = c(-.94,-.94,-.94), yend = c(-.06,-.06,-.06))
  to_phenos <- data.frame(x = c(-1,-1,.06), xend = c(-.06,1,1),
                          y = c(-1.06,-1.06,-1.06), yend = c(-1.06,-1.06,-1.06))
  
  paths <- data.frame(variables = letters[1:8],
                      x = c(-.5, .05, .5, -.55, .55, 0, .05, -.5),
                      y = c(-.55, -.55, -.55, -1.1, -1.1, -1.5, .3, .3))
  
  if(i==1){
    to_survive <- to_survive[2,]
    paths <- paths[paths$variables != "h",]
  }
  
  if(i==2){
    to_survive <- to_survive[1,]
    paths <- paths[paths$variables != "g",]
  }
  
  mods[[i]] <- data.frame(variable = c("Survive", "Size", "RGR", "Age","SeedSize"),
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




# plot size, rgr, size-rgr
plot_grid(plots[[1]], plots[[2]], plots[[3]], 
          mods[[1]], mods[[2]], mods[[3]],
          path_coefs[[1]], path_coefs[[2]], path_coefs[[3]],ncol = 3, labels = "AUTO")

ggsave(snakemake@output[[1]], device = "svg", width = 12, height = 8)

alpha <- fit$draws("alpha_survive", format = "df")$alpha_survive
beta <- fit$draws("beta_rgr_survive", format = "df")$beta_rgr_survive

data.frame(alpha, beta) %>% 
  pivot_longer(1:2, names_to = "coefficient", values_to = "estimate") %>% 
  group_by(coefficient) %>% 
  summarise(mean = mean(estimate),
            "97.5%" = quantile(estimate, .975),
            "75%" = quantile(estimate, .75),
            "25%" = quantile(estimate, .25),
            "2.5%" = quantile(estimate, .025)) #%>% 
  #write.csv(snakemake@output[['coef']], row.names = F, quote = F)

print("All plots made")

on.exit(sink())




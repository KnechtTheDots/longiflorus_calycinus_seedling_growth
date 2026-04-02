# the first three lines create the log file for the snakemake
# workflow
log <- file(snakemake@log[['log']], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(tidyr)
library(cmdstanr)

delta <- function(x){
  (x$f2 - (.5*x$f1 + .25*(x$lon + x$cal)))/x$e_f2
}

traits <- c("Survive", "Day_4", "Day_17", "RGR", "Height")
for(i in 1:5){
  d <- read.csv(snakemake@input[[i]])
  df <- d
  d$delta <- delta(d)
  
  # write the results of the epistasis test to a csv
  d <- d %>% 
    pivot_longer(1:ncol(d), names_to = "line", values_to = "value") %>% 
    group_by(line) %>% 
    summarise(quant.025 = quantile(value, .025),
              quant.25 = quantile(value, .25),
              mean = mean(value),
              quant.75 = quantile(value, .75),
              quant.975 = quantile(value, .975),
              prop_less_0 = sum(value < 0)/n()) %>%
    mutate(Trait = traits[i])
  
  # write the results of the pairwise differences between lines to a file
  dd <- data.frame("F2-F1" = df$f2 - df$f1,
             "F2-LON" = df$f2 - df$lon,
             "F2-CAL" = df$f2 - df$cal,
             "F1-LON" = df$f1 - df$lon,
             "F1-CAL" = df$f1 - df$cal,
             "LON-CAL" = df$lon - df$cal) %>% 
    pivot_longer(1:6, names_to = "comp", values_to = "value") %>% 
    group_by(comp) %>% 
    mutate(comp = case_when(comp == "F2.F1" ~ "F2-F1",
                            comp == "F2.LON" ~ "F2-LON",
                            comp == "F2.CAL" ~ "F2-CAL",
                            comp == "F1.LON" ~ "F1-LON",
                            comp == "F1.CAL" ~ "F1-CAL",
                            comp == "LON.CAL" ~ "LON-CAL"),
           comp = factor(comp, levels = c("F2-F1","F2-LON","F2-CAL",
                                          "F1-LON","F1-CAL","LON-CAL"))) %>% 
    summarise(quant.025 = quantile(value, .025),
              quant.25 = quantile(value, .25),
              mean = mean(value),
              quant.75 = quantile(value, .75),
              quant.975 = quantile(value, .975),
              prop_less_0 = round(mean(value < 0), 2)) %>% 
    mutate(trait = traits[i]) %>% 
    select(Trait = trait, Comparison = comp, quant.025, quant.25, mean, quant.75, quant.975, prop_less_0)
  
  if(i==1){
    d_traits <- d
    d_diff <- dd
  }else{
    d_traits <- rbind(d_traits,d)
    d_diff <- rbind(d_diff,dd)
  }
}
write.csv(d_traits, snakemake@output[['traits']], row.names = F, quote = F)
write.csv(d_diff, snakemake@output[['diffs']], row.names = F, quote = F)

# get the posteriors for the coefficients for the survival model and write to file
mod_size <- readRDS(snakemake@input[['mod_size']])
mod_rgr <- readRDS(snakemake@input[['mod_rgr']])
mod_rgr_size <- readRDS(snakemake@input[['mod_rgr_size']])

size <- data.frame(value = mod_size$draws("beta_size_survive", format = "df")$beta_size_survive, model = "Size Only", coefficient = "beta_size_survive")

rgr <- data.frame(value = mod_rgr$draws("beta_rgr_survive", format = "df")$beta_rgr_survive,
                  model = "RGR only", coefficient = "beta_rgr_survive")

rgr_comb <- data.frame(value = mod_rgr_size$draws("beta_rgr_survive", format = "df")$beta_rgr_survive, model = "Combined", coefficient = "beta_rgr_survive")

size_comb <- data.frame(value = mod_rgr_size$draws("beta_size_survive", format = "df")$beta_size_survive, model = "Combined", coefficient = "beta_size_survive")

d <- rbind(size, rgr, rgr_comb, size_comb)

d %>% 
  group_by(model, coefficient) %>% 
  summarise(quant.025 = quantile(value, .025),
            quant.25 = quantile(value, .25),
            mean = mean(value),
            quant.75 = quantile(value, .75),
            quant.975 = quantile(value, .975),
            prop_greater_0 = sum(value > 0)/n()) %>% 
  write.csv(snakemake@output[['betas']], row.names = F, quote = F)



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
            prop_less_0 = round(mean(diff < 0), 2))


# get the posteriors for the correlations
# function for calculating the path coefficients
get_cors <- function(r_rgr_size, r_age_size, r_seed_size, d, e, f){
  c <- ((r_seed_size - f*r_rgr_size) - (e - d*f)*(r_age_size - d*r_rgr_size)/(1-d^2))/
    (1 - f^2 - (e-d*f)^2/(1-d^2))
  b <- (r_age_size - d*r_rgr_size)/(1-d^2) - c*(e-d*f)/(1-d^2)
  a <- r_rgr_size - b*d - c*f
  return(data.frame(a, b, c, d, e, f))
}

fits <- list(mod_size, mod_rgr, mod_rgr_size)
models <- c("Size Only", "RGR Only", "Combined")

# read in data that goes into model
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


for(i in 1:3){
  
  fit <- fits[[i]]
  
  d <- fit$draws("R_phenos", format = "df")$`R_phenos[2,1]`
  e <- fit$draws("R_phenos", format = "df")$`R_phenos[2,3]`
  f <- fit$draws("R_phenos", format = "df")$`R_phenos[1,3]`
  
  rgr <- fit$draws("rgr", format = "df")[,1:nrow(short)]
  size <- fit$draws("final_rep", format = "df")[,1:nrow(short)]
  age <- dat$age_max
  seed <- dat$seed
  
  w <- fit$draws("w", format = "df")[,1:nrow(short)]
  sigma_w <- fit$draws("sigma_w", format = "df")$sigma_w
  
  r_rgr_size <- r_age_size <- r_seed_size <- r_rgr_surv <- r_size_surv <- c()
  for(j in 1:nrow(rgr)){
    r_rgr_size[j] <- cor(as.numeric(rgr[j,]), as.numeric(size[j,]))
    r_age_size[j] <- cor(age, as.numeric(size[j,]))
    r_seed_size[j] <- cor(seed, as.numeric(size[j,]))
    r_rgr_surv[j] <- cor(as.numeric(rgr[j,]), as.integer(w[j,]))
    r_size_surv[j] <- cor(as.numeric(size[j,]), as.integer(w[j,]))
  }
  
  cors <- get_cors(r_rgr_size, r_age_size, r_seed_size, d, e, f)
  
  if(i==1){
    cors$g <- r_size_surv
    cors$h <- 0
    cors$rgr_total <- cors$a * cors$g
  }
  
  if(i==2){
    cors$g <- 0
    cors$h <- r_rgr_surv
    cors$rgr_total <- cors$h
  }
  
  if(i==3){
    cors$g <- (r_size_surv - r_rgr_size*r_rgr_surv)/(1 - r_rgr_size^2)
    
    cors$h <- r_rgr_surv - cors$g*r_rgr_size
    
    cors$rgr_total <- cors$h + cors$a * cors$g
  }
  
  cors$s_rgr <- cors$rgr_total*sigma_w
  cors$s_age <- (cors$b * cors$g)*sigma_w
  cors$s_seedsize <- (cors$c * cors$g)*sigma_w
  cors$s_size <- cors$g * sigma_w
  
  cors <- cors %>% 
    pivot_longer(1:ncol(cors), names_to = "path", values_to = "value") %>% 
    group_by(path) %>% 
    summarise(quant.025 = quantile(value, .025),
              quant.25 = quantile(value, .25),
              mean = mean(value),
              quant.75 = quantile(value, .75),
              quant.975 = quantile(value, .975),
              prop_greater_0 = sum(value > 0)/n()) %>% 
    mutate(Model = models[i])
  
  if(i==1){
    cors_summary <- cors
  }else{
    cors_summary <- rbind(cors_summary, cors)
  }
    
}

write.csv(cors_summary, snakemake@output[['corrs']], quote = F, row.names = F)


on.exit(sink())
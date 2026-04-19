# the first three lines create the log file for the snakemake
# workflow
log <- file(snakemake@log[['log']], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(tidyr)

d <- read.csv(snakemake@input[['traits']])

d %>% 
  filter(line != "e_f2") %>% 
  group_by(Trait) %>% 
  mutate(line = factor(line, levels = c("cal", "lon", "f1", "f2", "delta"))) %>% 
  arrange(line, .by_group = T) %>% 
  mutate(line = case_when(line=="cal"~"CAL",
                          line=="lon"~"LON",
                          line=="f1"~"$F_1$",
                          line=="f2"~"$F_2$",
                          line=="delta"~"$\\Delta/E(F_{2})$"),
         Trait = gsub("_", " ", Trait)) %>% 
  write.csv(snakemake@output[['traits']], quote = F, row.names = F)

d <- read.csv(snakemake@input[['betas']])
d %>% 
  filter(model == "Combined") %>% 
  select(-model) %>% 
  mutate(coefficient = case_when(coefficient=="beta_rgr_survive"~"$\\beta_{rgr}$",
                                 coefficient=="beta_size_survive"~"$\\beta_{size}$")) %>% 
  write.csv(snakemake@output[['betas']], quote = F, row.names = F)

d <- read.csv(snakemake@input[['diffs']])

d %>% 
  mutate(Trait = gsub("_", " ", Trait),
         Comparison = gsub("F2", "$F_2$", Comparison),
         Comparison = gsub("F1", "$F_1$", Comparison)) %>% 
  write.csv(snakemake@output[['diffs']], quote = F, row.names = F)

d <- read.csv(snakemake@input[['corrs']])
d <- read.csv("results/posteriors/corrs.csv")
d %>% 
  filter(Model == "Combined",
         path != "rgr_total") %>% 
  mutate(path = case_when(path=="s_age"~"$\\beta_{age}$",
                          path=="s_rgr"~"$\\beta_{rgr}$",
                          path=="s_seedsize"~"$\\beta_{seed}$",
                          path=="s_size"~"$\\beta_{size}$",
                          T ~ path)) %>% 
  write.csv(snakemake@output[['corrs']], quote = F, row.names = F)


on.exit(sink())
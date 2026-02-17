# the first three lines create the log file for the snakemake
# workflow
log <- file(snakemake@log[['log']], open = "wt")
sink(log, type = "output")
sink(log, type = "message")


# define functions to do the bootstrapping and delta calculations
b_boot <- function(n = 1e4, trait){
  weights <- gtools::rdirichlet(n, rep(1, length(trait)))
  return(apply(weights, 1, function(x) sum(x*trait)))
}



boot_df_make <- function(f2, f1, lon, cal, trait){
  f2 <- b_boot(trait = f2[,trait])
  f1 <- b_boot(trait = f1[,trait])
  lon <- b_boot(trait = lon[,trait])
  cal <- b_boot(trait = cal[,trait])
  e_f2 <- .5*f1 + .25*(lon + cal)
  return(data.frame(f2, f1, lon, cal, e_f2))
}




traits <- read.csv(snakemake@input[['traits']])
traits$rgr <- (log(traits$day_17 - log(traits$day_4)))/(17 - 4)

# separate into data_frames by line

f2 <- traits[traits$line=="F2",]
f1 <- traits[traits$line=="F1",]
lon <- traits[traits$line=="LON",]
cal <- traits[traits$line=="CAL",]

# calculate survival probabilities
write.csv(boot_df_make(f2, f1, lon, cal, trait = "survive"),
          file = snakemake@output[['survive']],
          row.names = F, quote = F)

# calculate the day_4 means
write.csv(boot_df_make(f2, f1, lon, cal, trait = "day_4"),
          file = snakemake@output[['day_4']],
          row.names = F, quote = F)

# calculate day_17 means
write.csv(boot_df_make(f2, f1, lon, cal, trait = "day_17"),
          file = snakemake@output[['day_17']],
          row.names = F, quote = F)

# calculate rgr means
write.csv(boot_df_make(f2, f1, lon, cal, trait = "rgr"),
          snakemake@output[['rgr']],
          row.names = F, quote = F)

# calculate height means
write.csv(boot_df_make(f2[!is.na(f2$height),],
                       f1[!is.na(f1$height),],
                       lon[!is.na(lon$height),],
                       cal[!is.na(cal$height),],
                       trait = "height"),
          snakemake@output[['height']],
          row.names = F, quote = F)


on.exit(sink())


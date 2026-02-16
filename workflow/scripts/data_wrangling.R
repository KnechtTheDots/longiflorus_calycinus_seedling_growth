# the first three lines create the log file for the snakemake
# workflow
log <- file(snakemake@log[['log']], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(tidyr)

print("Packages loaded")

# read in the data with the day 4 seedling areas
d4 <- read.csv(snakemake@input[['day_4']])
# read in the data with the day 17 seedling areas
d17 <- read.csv(snakemake@input[['day_17']])
# read in the information that identifies the plants
germ_id <- read.csv(snakemake@input[['ids']])
# read in the height data
height <- read.csv(snakemake@input[['height']])

print("Data read in")

height <- height %>% 
  # the plant ids have periods instead of underscores so change from
  # periods to underscores
  mutate(plant = gsub("\\.", "_", plant)) %>% 
  # plant id is "six_pack" in the other data sets so rename it
  select(six_pack = plant, height) %>% 
  drop_na(six_pack)

# join the areas and height into a single data set
growth <- left_join(d4, germ_id, by = join_by(tray_id)) %>% 
  left_join(d17, by = join_by(six_pack)) %>% 
  drop_na(six_pack, day_4)

# some labels got either mistranscribed or mislabeled so there are duplicate individuals,
# I am dropping them.

dup_ids <- growth$six_pack[duplicated(growth$six_pack)]

growth <- growth %>% 
  filter(!(six_pack %in% dup_ids))

# repeat for height

dup_ids <- height$six_pack[duplicated(height$six_pack)]

height <- height %>% 
  filter(!(six_pack %in% dup_ids))

d <- left_join(growth, height, by = join_by(six_pack))


# filter so that only individuals that germinated on day 3 remain

d <- d %>% 
  filter(germ_day==3)



# a couple plants have missing data for size at day 17 but have data for height so they
# aren't dead, so I am removing them

# and by making height numeric, the dead individuals will be represented as NAs

d <- d %>% 
  drop_na(day_17) %>% 
  mutate(height = as.numeric(height),
         survive = ifelse(is.na(height), 0, 1)) %>% 
         select(line, day_4, day_17, height, survive)

write.csv(d, snakemake@output[['traits']], row.names = F, quote = F)

print("Traits data written to file")

# read in the drought trial data
seed_size <- read.csv(snakemake@input[['drought']])

# create the 'long' data set where each row is a measurement of a plant
# at a specific timepoint
long <- seed_size %>% 
  drop_na(germ_day) %>% 
  pivot_longer(day_3:day_17, names_to = "day", values_to = "seedling_area") %>% 
  drop_na(seedling_area) %>% 
  mutate(day = as.numeric(gsub("day_","", day)),
         age = day - germ_day) %>% 
  # since it is hared to measure them on the germ day, they are still unfolding,
  # I am removing measurements from that day
  filter(age > 0) %>% 
  drop_na() %>%
  # the last day for 1_22 must be measured wrong, need to go back and
  # re-measure it, but for now remove it.
  filter(id != "1_22") %>% 
  mutate(id = as.integer(factor(id)),
         id = factor(id, levels = 1:length(unique(id))))

write.csv(long, snakemake@output[['long']], row.names = F, quote = F)

# create the 'short' dataset where each row is the final size and 
# survival of a plant during the drought trial.
short <- long %>% 
  group_by(id) %>% 
  filter(age == max(age)) %>% 
  # did they survive till day 4? (1,0)
  mutate(survive = death_day - 3) %>% 
  arrange(id) %>% 
  select(id, seed_size, seedling_area, survive, age)

write.csv(short, snakemake@output[['short']], row.names = F, quote = F)

print("Drought trial data written to file")


on.exit(sink())
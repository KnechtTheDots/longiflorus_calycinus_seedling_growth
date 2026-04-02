The csv's in this directory contain the results for two experiments. The first, described as epistasis test below contains
the trait measurements for survival, size on day 4, size on day 17, height on day 122, and relative growth rate (n.b. survival and rgr are
calculated as part of the analysis pipeline and aren't contained in one of the source files below). The second, described as the survival/drought trial
contains information on traits and how they relate to fitness in a terminal drought experiment.

## Files used for the test for epistasis are:

**day_4_areas.csv** with the following columns  
*tray_id*: this is the id of each individual plant  
*day_4*: area (in mm^2) of the seedling on day 4 of the experiment  

**day_17_areas.csv** with the following columns:  
*six_pack*: this is the id of each individual plant, at this point they have been transplanted from trays into six packs  
*day_17*: area (in mm^2) of the seedling on day 17 of the experiment  

**heights.csv** with the following columns:  
*plant*: this is the id of each individual plant--matching the six pack id in the day_17 csv  
*height*: the height (to the nearest 0.5 cm) of each plant  

**germ_and_id.csv** with the following columns:  
*tray_id*: the tray id of each plant  
*six_pack*: the six_pack id of each plant (i.e. links the ids pre- and post transplant)  
*germ_day*: the day the plant germinated (from the start of the experiment)  
*cross*: identification of the parents format=motherxfather  
*line*: which line the individual belongs to (i.e. lon, cal, f1, or f2)  
*notes*: misc. notes  

## Files used for the drought trial  

**seed_growth_vs_seed_size.csv** with the following columns:  
*id*: the id of each seed  
*seed_size*: the size (area in mm^2) of each seed  
*germ_day*: the day (time since beginning the experiment) the seed germinated  
*day_3-17*: size of seedling (area in mm^2) on each day measured during the experiment  
*death_day*: day of death of each seedling in the drought trial (time since the beginning of the drought trial)  
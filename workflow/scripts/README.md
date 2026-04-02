This directory contains all the scripts necessary for the analysis.  

## The **stan** sub-directory contains 3 stan programs that fit models relating growth rate, age, and seed size to fitness.  

*stan/drought_survival_rgr.stan* models fitness as a function of growth rate.  
*stan/drought_survival_size.stan* models fitness as a function of seedling size.  
*stan/drought_survival_rgr_size.stan* models fitness as a function of both size and growth rate.  
*compile_complete.txt* is a dummy file that is generated as part of the workflow (See the "compile_stan" rule in *Snakefile* for more details)  

## R scripts

*data_wrangling.R* takes the raw data files and cleans them and creates output files for analysis.  
*trait_boots.R* creates files containing the bayesian bootstraps for the epistasis test  
*epi_plots.R* creates plots for the results for the epistasis test as well as comparing crosses and plotting the raw trait values    
*drought_survival.R* runs the stan models and creates RDS files of the results  
*ppc_plots.R* conducts posterior predictive and diagnostic checks to validate the stan models and generates a plot  
*mod_compare.R* does model comparison of the fitted models  
*survival_plots.R* creates plots to show the results of the stan models  
*get_posteriors.R* creates csv's describing the posteriors for all the quantities of interest
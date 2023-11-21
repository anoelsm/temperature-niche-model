#### Smith and Barton (2023) Metacommunity model 
### Version: F76F9sFK8xZa
### scp -r ~/Documents/PhD_Work/Research/Niche_Model/Code/{diversity_function-update.R,tniche_analysis-update.R} ans132@kaweah.ucsd.edu:documents/metacomm_model/
#### scp -r ans132@kaweah.ucsd.edu:documents/metacomm_model/csv_output/{converged_values_ensemble.csv, ~/Documents/PhD_Work/Research/Niche_Model/Data/
#### scp -r ans132@kaweah.ucsd.edu:documents/metacomm_model/ncdf/reviews/review_analysis-ensemble.Rdata ~/Documents/PhD_Work/Research/Niche_Model/Data/
#### scp -r ans132@kaweah.ucsd.edu:documents/metacomm_model/ncdf/final/* ~/Documents/PhD_Work/Research/Niche_Model/Data/nc_files/

#### Packages ----
packages <- c(
  "ncdf4", "tidyr", "class", "fields", "R.matlab", "vegan",
  "gdata", "propagate", "nlstools", "stringi", "plyr", "sm",
  "dplyr", "zoo", "gam", "tidyr", "reshape", "stringr", "rlang",
  "tidyverse", "abind", "rbin", "plyr", "io", "rlang", "purrr")

funlist <-  lapply(packages, function(x) {
  if (x %in% rownames(installed.packages())) {
    library(x, character.only = T)
  }else{
    install.packages(x, character.only = T); library(x, character.only = T)
  }
})

# Personal preference to max out console output
options(max.print = 100)

### Functions ---
source("documents/metacomm_model/tniche_analysis-update.R")
source("documents/metacomm_model/diversity_function-update.R")

#### Folders ----
files_home <- c("documents/metacomm_model/ncdf/reviews/")
files_home <- c("documents/metacomm_model/ncdf/")
nc_path <- c("~/Documents/PhD_Work/Research/Niche_Model/Data/nc_files/")
#### Define variables -----
sp_sep <- 1 #Number of degrees between optimal temperatures
z_vals <- seq(-4, 40, by = sp_sep) #Create the range of possible z values
num_species <- length(z_vals) #Get max number of species
spacing <- 1 #Number of degree latitude between boxes
full_latitude <- seq(-79, 79, by = spacing) # Discretation of full_latitudes
latitude <- seq(-79, 79 + spacing, by = spacing) #Add one to round it out for ranges
bx <- length(full_latitude) #Max number of boxes
boxes_to_plot <- which(full_latitude %in% c(10, 35, 60, 75)) #Pull out 4 illustrative boxes
k_h <- c(c(1E0, 1E1, 1E2, 1E3)) * (60 * 60 * 24) #Dispersal magnitudes
lats_plot <- seq(-75, 75, length.out = 7) #Latitudes to plot of y axis
years <- 50 #How many years the model runs for; should up to 100 on final runs
n_day <- 365 * years #Number of days for model run purposes
t_step <- 3/24 #Number of hours to step forward in model
iter <- n_day / t_step #Max model iterations
oneyr_iter <- 365 / t_step #How many model iterations in one year
time_lastyr <- ((oneyr_iter*5) - oneyr_iter):(oneyr_iter*5) #Model iterations for the last year
oneyr_hours <- seq(1,366-t_step, by = t_step)
lat_t_mean <- matrix(
  nrow = length(ttime),
  ncol = bx) #Initiate a matrix with rows = temperature time and col = number of boxes
# For each latitude band, average the temperatures across all longitudes at each time point
for (i in seq_along(full_latitude)){
  latband <- which(
    tlat >= latitude[i] & tlat <= latitude[i+1])
  lat_t_mean[,i] <- apply(
    tsst[,latband,], 3, function(x){
      mean(x,na.rm=T)})
} 
# Repeat the seasonal cycle for the number of years
trep <- do.call("rbind", rep(list(lat_t_mean), years))
# Discretize the time steps (1 day) to the model time steps (3 hours)
time_vec <- seq(1, n_day, len = iter) #time_vec vector
mlist <- apply(
  trep, 2, function(x) {approx(x, xout = time_vec)})
# Add all to one data frame
temp_mat <- do.call(
  "cbind",lapply(mlist, function(c){c$y}))

#### Define fundamental niches -----
fund_vals <- lapply(1:num_species, function(x) {
  x.vals <- seq(-40, 50, by = 0.1)
  y.vals <- 0.81 * exp(0.0631 * x.vals) *
    (1 - ((x.vals - z_vals[x]) / (10 / 2))^2)
  y.vals[y.vals < 0] <- 0
  return(list(
    vals = data.frame(
      topt = x.vals[which.max(y.vals)],
      width = diff(range(x.vals[y.vals > 0]))),
    df = data.frame(
      x = x.vals,
      y = y.vals,
      sp = x)))
})
### Pull out the widths and the topts seperately
fund_widths <- unlist(lapply(
  fund_vals, function(x) x$vals$width), use.names = F)
fund_topt <- unlist(lapply(
  fund_vals, function(x) x$vals$topt), use.names = F)
fund_curves <- do.call("rbind", lapply(
  fund_vals, function(x) x$df))

# file, subtitle, wdir, nc_path
nc_path <- "./documents/metacomm_model/ncdf/final/"
wdir <- "./documents/metacomm_model/csv_output/"
all_files <- list_files(nc_path)
all_analyses <- lapply(all_files, tniche_analysis, subtitle = "all_redo_final_2", wdir = wdir, nc_path = nc_path)
# dir.create("./documents/metacomm_model/ncdf/final/", recursive=TRUE)
# lapply(all_files, function(f){
#   file.rename(from = paste0("./documents/metacomm_model/ncdf/final", f),  to = paste0("./documents/metacomm_model/ncdf/final/", f))
# })

save(all_analyses, file = "./documents/metacomm_model/csv_output/all_diversity_list-publication-4.Rdata")

#scp -r ans132@kaweah.ucsd.edu:documents/metacomm_model/csv_output/{converged_values_all_redo_final_2.csv,all_diversity_list-publication-4.Rdata}  ~/Documents/PhD_Work/Research/Niche_Model/Data/

#scp ans132@kaweah.ucsd.edu:documents/metacomm_model/ncdf/T1_D0_E0_exp3-publication-2.nc ~/Documents/PhD_Work/Research/Niche_Model/Data/
file <- all_files[4]
range(all_analyses[[6]]$div$ids)
  
  
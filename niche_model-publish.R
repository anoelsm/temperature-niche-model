### Smith and Barton (2023) Metacommunity model 
### Version: F76F9sFK8xZa
### scp -r ~/Documents/PhD_Work/Research/Niche_Model/Code/tniche_analysis-update.R ans132@kaweah.ucsd.edu:documents/metacomm_model/
### scp -r ans132@kaweah.ucsd.edu:documents/metacomm_model/csv_output/all_diversity_list-publication-mini.Rdata ~/Documents/PhD_Work/Research/Niche_Model/Code/
#scp -r ans132@kaweah.ucsd.edu:documents/metacomm_model/ncdf/T0_D0_E0_exp1-review_ensemble.nc ~/Documents/
packages <- c("tidyr", "stringi", "plyr", "ncdf4",
"dplyr", "zoo", "stringr", "scales", "reshape")

funlist <- lapply(packages, function(c){
  if(!require(c, character.only = TRUE)){
    install.packages(c)
    require(c, character.only = TRUE)
  }else{
    require(c, character.only = TRUE)
  }
})
options(max.print = 100) #99999

##### Folders --------
Rds.outputfiles <- c("~/documents/metacomm_model/rds.output/")
nc_path <- "./documents/metacomm_model/ncdf/"

#Rds.outputfiles <- c("~/documents/rds.output/final/rerun/")
#### Load in climatological temperature data ----
load(paste0("./documents/metacomm_model/temps.Rdata"))

### Initial parameters ----
subtitle <- "-review_ensemble"

sp_sep <- 1
z_vals <- seq(-4, 40, by = sp_sep)
num_species <- length(z_vals) #species
t_step <- 3 / 24 #day
k <- 1 * (1 / 1E-3) * (1 / 1E3) #mmolPm^-3 half-saturation constant
d <- (1E-5) * 60 * 60 * 24 #vertical flux of nutrients from the deep; day^-1
R_0 <- (0.8) * (1 / 1E3) * (1 / 1E-3) #Concentration of nutrients in the deep being advected up; chemostatic mmolP m^-3
years <- 50 #How many years the model runs for; should up to 100 on final runs
n_day <- 365 * years #Number of days for model run purposes
iter <- n_day / t_step #Model iterations
P_alpha <- 1
spacing <- 1
full_latitude <- seq(-79, 79, by = spacing) ### Discretation of full_latitudes
latitude <- seq(-79, 79 + spacing, by = spacing) #add one to round it out
bx <- length(full_latitude)

boxes_to_plot <- which(full_latitude %in% c(10, 35, 60, 75))
k_h <- c(c(1E0, 1E1, 1E2, 1E3)) * (60 * 60 * 24)

time_vec <- seq(1, n_day, len = iter) #Time vector

oneyr_iter <- 365 / t_step
oneyear <- oneyr_iter
lastyear <- (iter - oneyr_iter):iter
lastyears <- (iter - oneyr_iter * 5):iter
yr_45 <- (365*45)/t_step
####### Take the average of the climatological sst for each latband 
lat_t_mean <- matrix(nrow = length(time_vec), ncol = bx) #For each box set up an empty vector of length time

#From the data file, take the mean temperature between each fulllatitude band
for (i in 1:length(full_latitude)){
  latband <- which(tlat >= latitude[i] & tlat <= latitude[i+1])
  lat_t_mean[,i] <- apply(tsst[,latband,], 3, function(x){mean(x, na.rm = T)})
}

#Repeat the yearly temp cycle for the number of years set
trep <- do.call("rbind", rep(list(lat_t_mean), years))
#Discretize the temperature values to the time vector
mlist <- apply(trep, 2, function(x){approx(x, xout = time_vec)})
temp_in <- do.call("cbind", lapply(mlist, function(c){c$y}))
temp_0 <- matrix(apply(temp_in, 2, mean), nrow = nrow(temp_in), ncol = ncol(temp_in), byrow = TRUE) 

#### Setup fundamental niches ----
mu_1 <- array(data=NaN, dim = c(length(time_vec), num_species, bx)) 
mu_0 <- mu_1 

z <- matrix(rep(z_vals, bx), bx, num_species, byrow = TRUE)
w <- matrix(10, bx, num_species, byrow = TRUE) 

for (j in 1:num_species){
  mu_1[,j,] <- 0.81*exp(0.0631*temp_in)*(1-((temp_in-z[,j])/(w[,j]/2))^2)
  mu_0[,j,] <- 0.81*exp(0.0631*temp_0)*(1-((temp_0-z[,j])/(w[,j]/2))^2)
}
mu_1[mu_1 < 0] <- 0 
mu_0[mu_0 < 0] <- 0 

#### Setting up dispersal ----
######Create an empty matrix that is bx by bx dimensions
taumatrix_0 <- matrix(0, nrow = bx, ncol = bx)
######Name the rows and columns amu.ter the box numbers
colnames(taumatrix_0) <- full_latitude
rownames(taumatrix_0) <- full_latitude
######Take the difference between the rows and columns for the difference factor
diff_matrix <- abs(outer(as.numeric(colnames(taumatrix_0)), as.numeric(rownames(taumatrix_0)), "-"))
colnames(diff_matrix) <- full_latitude
rownames(diff_matrix) <- full_latitude
####Give realistic differences (1 degree fulllatitude= 111km convert to m)
diff_matrix <- diff_matrix*(111*1000)

###### Run Model-----

## Set up netCDF

nc_species <- ncdim_def("species", "ind", 1:num_species)
nc_lat <- ncdim_def("latitude", "degrees", full_latitude)
nc_time <- ncdim_def("time", "days", time_vec[yr_45:iter])

nc_mu <- ncvar_def("fund_niche", "per_day", list(nc_time,nc_species,nc_lat))
nc_p <- ncvar_def("biomass", "mmolP_per_m3", list(nc_time,nc_species,nc_lat))
nc_r <- ncvar_def("nutrients", "mmolP_per_m3", list(nc_time,nc_species,nc_lat))
nc_growth <- ncvar_def("net_growth", "mmolP_per_m3_per_day", list(nc_time,nc_species,nc_lat))
nc_decay <- ncvar_def("net_decay", "mmolP_per_m3_per_day", list(nc_time,nc_species,nc_lat))
nc_tausum <- ncvar_def("tau_in", "mmolP_per_m3_per_day", list(nc_time,nc_species,nc_lat))
nc_tausub <- ncvar_def("tau_out", "mmolP_per_m3_per_day", list(nc_time,nc_species,nc_lat))
nc_temp <- ncvar_def("temp", "degrees_C", list(nc_time,nc_lat))

P <- array(
  data = NA,
  dim = c(iter, num_species, bx),
  dimnames = list(
    paste0("t_", 1:iter),
    paste0("species_", 1:num_species),
    paste0("lat_", signif(full_latitude, 4))))
R <- P
tausum <- P
tausub <- P
growth <- P
decay <- P

P[1,,i] <- 1E-3
R[1,,] <- 1E-3



for (wm in 1:3){
  m <- c(1/20, 1/10, 1/5)[wm]
  print(paste("I'm tired of this grandpa ~ ", wm))
  for (e in 1:2){
    #Model Conditions
    variable <- c(0, 1, 0, 1)[e]
    interaction <- c(0, 0, 1, 1)[e]
    if (variable == 0){
      temp <- temp_0
      mu_t <- mu_0
    }else{
      temp <- temp_in
      mu_t <- mu_1
    }
    t_mortality <- (0.81 * exp(0.0631 * temp))
    
    #mu_t[(iter/2):iter,20,c(1:39,41:bx)] <- 0
   
    for (kval in seq_along(k_h)){
      print(paste("I can fix that ~ ", kval))

      taumatrix_2 <- k_h[kval]/(diff_matrix^2)
      taumatrix_2[which(!is.finite(taumatrix_2))] <- 0

      ##If there is no interaction then the taumatrix is just zeros
      if (interaction == 1){
        taumatrix_2 <- taumatrix_2
      }else {
        taumatrix_2 <- taumatrix_0
        kval <- 0
      }

      ##Create the filename within this loop
      output_filename <- paste0("T", variable, "_D", interaction, "_E", kval, "_exp", wm, subtitle)
      nc_fname <- paste(nc_path, output_filename, "-2.nc", sep="")
      nc_out <- nc_create(
        nc_fname,
        list(nc_mu, nc_p, nc_r, nc_growth, nc_decay, nc_tausum, nc_tausub, nc_temp),
        force_v4 = TRUE)
      
      ncatt_put(nc_out, 0, "title", paste("Last 5 years of metacommunity model run experiment -", output_filename))
      ncatt_put(nc_out, 0, "institution","Scripps Institution of Oceanography (UCSD)")
      ncatt_put(nc_out, 0, "authors", "Alaina N. Smith & Andrew D. Barton")
      history <- paste("A.N. Smith", date(), sep=", ")
      ncatt_put(nc_out,0, "history", history)
      
      ncvar_put(nc_out, nc_mu, mu_t[yr_45:iter,,])
      ncvar_put(nc_out, nc_temp, temp[yr_45:iter,])
      
      ####### Initial Values ----
     
      ###### Run the model ----
      for (i in 2:iter){
        if (i == 2){
          print("What does D-I-G spell?")
        }
        
        P2 <- array(apply(P[i-1,,], 1, function(x){rep(x, each = bx)}), dim = c(bx, bx, num_species)) 
        
        tausum[i-1,,] <- t(apply(P2, 3, function(x){
          apply(x*taumatrix_2, 1, function(y){
            sum(y, na.rm = T)})}))
        tausub[i-1,,] <- t(apply(P2, 3, function(x){
          apply(x*taumatrix_2, 2, function(y){
            sum(y, na.rm = T)})}))
        
        growth[i-1,,] <- mu_t[i-1,,] * P[i-1,,] * (R[i-1,,]/(R[i-1,,]+k))
        const_mortality <- m*(P[i-1,,]^P_alpha)
        decay[i-1,,] <- sweep(const_mortality,2,t_mortality[i-1,],FUN='*')
        chemostat <- d*(R_0-R[i-1,,])
        
        P[i,,] <- P[i-1,,]+t_step*(growth[i-1,,] - decay[i-1,,] + tausum[i-1,,] - tausub[i-1,,])
        R[i,,] <- R[i-1,,]+t_step*sweep(chemostat, 2, apply(growth[i-1,,], 2, function(x){
          sum(x,na.rm=T)
            }), FUN = '-')
        
        if(i == oneyr_iter*30){
          #After 30 years, save the last output
          p_mid <- P[i,,]
          r_mid <- R[i,,]
          save(P_init = p_mid, R_init = r_mid,
            file = paste0(Rds.outputfiles, output_filename, "-year30.Rdata"))
          print(paste0("30 years: ", output_filename))
        }
      }
      
      ncvar_put(nc_out,nc_p,P[yr_45:iter,,])
      ncvar_put(nc_out,nc_r,R[yr_45:iter,,])
      ncvar_put(nc_out,nc_growth,growth[yr_45:iter,,])
      ncvar_put(nc_out,nc_decay,decay[yr_45:iter,,])
      ncvar_put(nc_out,nc_tausum,tausum[yr_45:iter,,])
      ncvar_put(nc_out,nc_tausub,tausub[yr_45:iter,,])
      
      nc_close(nc_out)
      print("It spells DIG!")

      if (interaction == 0){
        break
      }
    }
  }
}

#### Pau -----

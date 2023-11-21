
# Download files from github
### Set up ----
#scp ans132@kaweah.ucsd.edu:documents/metacomm_model/csv_output/ ~/Documents/PhD_Work/Research/Niche_Model/Data/
#### Download packages ----
packages <- c(
  "ncdf4", "tidyr", "RColorBrewer", "ggplot2", "class", "fields", "R.matlab",
  "caTools", "readxl", "gdata", "propagate", "nlstools", "stringi", "plyr",
  "XML", "biomod2", "grid", "pwr", "formattable", "dplyr", "zoo", "lattice",
  "gridExtra", "gam", "tidyr", "cowplot", "nls2", "reshape", "stringr",
  "tidyverse", "abind", "rbin", "patchwork", "ggpubr", "sm",  "io")

funlist <-  lapply(packages, function(x) {
  if (x %in% rownames(installed.packages())) {
    require(x, character.only = T)
  }else{
    install.packages(x, character.only = T); require(x, character.only = T)
  }
})

# Figure settings 
colors_P <- colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red"))
niche_colors <- c("salmon1", "plum3", "palevioletred3")
mort <- c(1 / 20, 1 / 10, 1 / 5)
names(mort) <- paste("m_", seq_along(mort), sep = "")
mort_labs <- do.call("expression", lapply(
  mort, function(c) {
    substitute(italic(gamma) == X, list(X = c))}))
labs_kh <- do.call("expression", lapply(0:3, function(r) {
  if (r == 0){
    "10^0\nLow dispersal"}
  else if (r == 3){
    "10^3\nHigh dispersal"
  }
  else{
    substitute(10^X, list(X = r))
  }
}))

# Personal preference to max out console output
options(max.print = 100)

#### Set up directories ----
files_home <- "~/Documents/PhD_Work/Research/Niche_Model/"
fig_dir <- paste0(files_home, "/Figures/Reviews/")
nc_files <- paste0(files_home, "Data/nc_files/")

#### Load in model parameters (number of species, latitudes, temperatures) -----
source_files <- paste0(files_home, "Code/", c("Smith&Barton-model_parameters-2023.R", "diversity_function-update.R", "tniche_analysis-update.R"))
lapply(source_files, source)

#### Option A: Run model analysis----
files_home <- "./documents/metacomm_model/"
nc_path <- paste0(files_home, "Data/nc_files/") #"./documents/metacomm_model/ncdf/"
all_files <- list_files(nc_path)
all_analyses <- lapply(all_files, tniche_analysis, subtitle = "all_redo_final_5",
                       wdir = paste0(files_home, "csv_output/"), nc_path = nc_path)
save(all_analyses, file = paste0(files_home, "rds.output/all_diversity_list-publication-7.Rdata"))

#### Option B: Load analysis ----
converged_vals <- read_csv(paste0(files_home,"Data/converged_values-2.csv"))
load("~/Documents/PhD_Work/Research/Niche_Model/Data/all_diversity_list-publication-2.Rdata")

### Fig. 1 - Model Schematic -----
#### Built in Adobe Illustrator
### Fig. 2 - Hypotheses ----
#### Built in Adobe Illustrator
### Fig. 3, 4, 7, 10 - E1-4 Model Output ----

# Pull out the illustrative experiments
df_id <- c("T0_D0_E0_exp1", "T0_D1_E3_exp1", "T1_D0_E0_exp1", "T1_D1_E3_exp1")
# Pull out all possible ids from diversity list
id_list <- do.call("rbind",(lapply(all_analyses,function(f){f$div$ids})))
# Pick a species to use for Fig. 11
sp_sink <- 24 
# Loop through each experiment and pull out the temperature and biomass timeseries for
# each box and each species that has converged in the model
for (e in 1:4){
  # open nc file
  nc_file_name <- list_files(nc_files)[grepl(df_id[e], list_files(nc_files))]
  nc <- nc_open(paste0(nc_files, nc_file_name))

  # identify which list element matches the variables
  exp_id <- which(apply(id_list, 1, function(r){all(str_split(df_id[e], "_")[[1]] %in% r)}))
  
  ##### For sink species niche width and source-sink dynamics -----
  if (e >= 2){
    # Determine which boxes one species (sp_sink) is present in
    box_pres <- which(unlist(lapply(all_analyses[[exp_id]]$div$species, function(r){
      any(sp_sink %in% r)})))
    # For each box, find out which column species is found in "total_alpha"
    column_sp <- unlist(lapply(all_analyses[[exp_id]]$div$species, function(r){
      which(sp_sink == r)}))
    # For each box, pull out the presense/absence time series for sp_sink,
    # and create a dataframe showing which box, the temperature in that box,
    # the biomass timeseries for sp_sink in that box,
    # wether sp_sink was present at each iteration or not,
    # and the id variables
    # then remove rows where sp_sink is not present
    converged_df_sp <- do.call("rbind",lapply(seq_along(box_pres), function(r){
      alpha_pres <- all_analyses[[exp_id]]$div$total_alpha[[
        box_pres[r]]][
          time_lastyr,column_sp[r]]
      df <- data.frame(
        box = box_pres[r],
        temp = as.vector(ncvar_get(nc, "temp")[time_lastyr,box_pres[r]]),
        bio =  ncvar_get(nc, "biomass")[time_lastyr, sp_sink, box_pres[r]],
        pres = c(alpha_pres),
        id = df_id[e]
      )
      return(df)
    })) %>%
      subset(., pres == TRUE)
    
    if (e == 4){
      # During the final experiment (E4 - combined) pull out all of the information for
      # sp_sink in each box and save it into a dataframe
      boxes <- c(1:bx)[unlist(lapply(all_analyses[[exp_id]]$div$species, function(y){
        sp_sink %in% y}))]
      p_sink_df <- data.frame(
        time_vec = rep(time_lastyr/oneyr_iter, bx),
        temp = as.vector(ncvar_get(nc, "temp")[time_lastyr,]),
        bio = matrix(ncvar_get(nc, "biomass")[time_lastyr,sp_sink,]),
        nutr = matrix(ncvar_get(nc, "nutrients")[time_lastyr,sp_sink,]),
        t_mu = matrix(ncvar_get(nc, "fund_niche")[time_lastyr,sp_sink,]),
        t_mort = (1/20)*(0.81 * exp(0.0631 * as.vector(ncvar_get(nc, "temp")[time_lastyr,]))),
        net_growth =  matrix(ncvar_get(nc, "fund_niche")[
            time_lastyr,sp_sink,]) - (1/20)*(0.81 * exp(0.0631 * as.vector(
              ncvar_get(nc, "temp")[time_lastyr,]))),
        decay = matrix(as.vector(ncvar_get(nc, "net_decay")[time_lastyr,sp_sink,])),
        growth = matrix(as.vector(ncvar_get(nc, "net_growth")[time_lastyr,sp_sink,])),
        trans = matrix(as.vector(ncvar_get(nc, "tau_in")[time_lastyr,sp_sink,]) - as.vector(ncvar_get(nc, "tau_out")[time_lastyr,sp_sink,])),
        box = rep(1:bx, each = length(time_lastyr))) 
    }
    
    if (e == 2){
      # After E2, save the dataframe of the converged species in each box as a copy df
      converged_df_sp2 <- converged_df_sp
    }else{
      # Attach the dataframe of the converged species in each box to the copy df
      converged_df_sp2 <- rbind(converged_df_sp2, converged_df_sp)
    }
  }
  ##### For intro figures -----
  #Pull the surviving species in each box
  pres_sp <- lapply(boxes_to_plot, function(r){
    unlist(all_analyses[[exp_id]]$div$species[[r]])})
  # Combined the time, temperature, and biomass for each surviving species
  # in each box into one dataframe
  full_df <- do.call("rbind",lapply(seq_along(boxes_to_plot), function(r){
    data.frame(
      time_vec = rep(oneyr_hours, 5),
      temp = matrix(ncvar_get(nc, "temp")[1:(oneyr_iter*5),boxes_to_plot[r]]),
      bio = matrix(ncvar_get(nc, "biomass")[1:(oneyr_iter*5),pres_sp[[r]],boxes_to_plot[r]]),
      sp = rep(pres_sp[[r]], each = (oneyr_iter*5)),
      box = rep(boxes_to_plot[[r]], length(oneyr_iter*5)))
  }))
  # Average the temperature and biomass for each species in each box from each day of the year over the last 5 years
  plot_df_mean <- full_df %>%
    mutate(time_day = floor(time_vec)) %>%
    group_by(sp,box,time_day) %>%
    summarize(
      temp_mean = mean(temp),
      bio_mean = mean(bio)
    )
  
  # Pull out the unique species
  sp_names <- unique(plot_df_mean$sp)[order(unique(plot_df_mean$sp))]
  # Pull out the fundamental optimal temperatures for each unique species
  sp_use_plot <- fund_topt[as.numeric(sp_names)][order(
    fund_topt[as.numeric(sp_names)])]
  names(sp_use_plot) <- sp_names[order(as.numeric(sp_use_plot))]
  
  # Add in column for linetype (lt) and experiment id (exp)
  plot_df_mean$lt <- "np"
  plot_df_mean$exp <- e
  
  
  if (e == 1) {
    # For the first (null) experiment (E1) pull out the species, box, and make linetype "p"
    sp_names_none <- plot_df_mean %>%
      dplyr::select(sp,box) %>%
      distinct() %>%
      mutate(lt = "p")
    # Save a copy dataframe
    plot_df2 <- plot_df_mean
  } else {
    # For each additional experiment, in each box, make the linetype for the species present
    # in E1 (sp_names_none) "p" and leave the rest "np" 
    for (bt in 1:4){
      plot_df_mean$lt[
        plot_df_mean$box == sp_names_none$box[bt] &
          plot_df_mean$sp == sp_names_none$sp[bt]] <- "p"
    }
    # Add the dataframe to the copy dataframe
    plot_df2 <- rbind(plot_df2, plot_df_mean)
  }
  
  # Close the nc file
  nc_close(nc)
  
}

##### Make figures -----
# Reverse the order of boxes so that the warmest is at the bottom
plot_df2$box <- factor(plot_df2$box, levels = rev(boxes_to_plot))
# Set up tag labels (a., b., c., etc.)
labels_temp <-paste0(letters[c(1,3,5,7)], ".")
names(labels_temp) <- rev(boxes_to_plot)
labels_bio <- paste0(letters[c(2,4,6,8)], ".")
names(labels_bio) <- rev(boxes_to_plot)

plot_df2$lt[plot_df2$exp == 1] <- "p"
plot_df2$sp <- factor(plot_df2$sp, levels = unique(plot_df2$sp)[order(unique(plot_df2$sp))])
sp_colors <- colors_P(length(unique(plot_df2$sp)))
sp_use_plot <- fund_topt[as.numeric(levels(plot_df2$sp))]
names(sp_use_plot) <- levels(plot_df2$sp)
names(sp_colors) <- levels(plot_df2$sp)

intro_label_temp <- data.frame(
  box = boxes_to_plot,
  label = paste0(full_latitude[boxes_to_plot], "°N")
)

plot_df2$lt[plot_df2$exp == 1] <- "np"

#plot_df2$colors <- factor(plot_df2$sp, levels = names(sp_colors), labels = sp_colors)
for (g in 1:4){
  plot_sub <- subset(plot_df2, exp == g)
  plot_a <- ggplot(plot_sub)+
    geom_line(aes(x = time_day, y = temp_mean)) +
    scale_y_continuous(name = bquote(Temperature~(degree * C))) +
    scale_x_continuous(
      name = bquote(Day~of~the~year),
      expand = c(0,0),
      breaks = c(1, 90, 180, 270, 360),
      labels = c(1, 90, 180, 270, 360)) +
    geom_text(data = intro_label_temp, aes(x = -180, y = 15, label = label),
              fontface = "bold", size = 4.25, angle = 90) +
    facet_wrap(box ~ ., ncol = 1,
               labeller = labeller(box = labels_temp)) +
    theme_bw() +
    ggtitle("Temperature") +
    coord_cartesian(ylim = c(-1.5,29), xlim = c(0,365), clip = "off") +
    theme(
      plot.margin = unit(c(-0.1,0.1,0,1), "lines"),
      panel.grid = element_blank(),
      text = element_text(color = "black", size = 10),
      panel.border = element_rect(color = "black", linewidth = 0.3),
      legend.title = element_text(size = 10),
      strip.background = element_blank(),
      aspect.ratio = 1,
      panel.spacing.y = unit(0, "cm"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9, color = "black"),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold", vjust = -5),
      strip.text = element_text(hjust = -0.03, vjust = 0, face = "bold", size = 10))
      
  plot_b <- ggplot(plot_sub) +
    geom_line(aes(
      x = time_day,
      y = log10(bio_mean),
      col = factor(sp),
      group = factor(sp),
      linetype = lt)) +
    scale_color_manual(
      name = bquote(italic(T[F]^opt) ~ (degree * C)),
      labels = signif(sp_use_plot,3),
      values = sp_colors,
      guide = guide_legend(reverse = TRUE, ncol = 1, keyheight=0.3)) +
    scale_linetype_manual(
      values = c("np" = "solid", "p" = "dotdash"),
      guide = "none") +
    scale_y_continuous(
      name = bquote(Biomass~(log[10]~mmol~P~m^-3))) +
    scale_x_continuous(
      name = bquote(Day~of~the~year),
      expand = c(0,0),
      breaks = c(1, 90, 180, 270, 360),
      labels = c(1, 90, 180, 270, 360)) +
    coord_cartesian(ylim = c(-6, 1.5), xlim = c(0, 365), expand = F) +
    facet_wrap(. ~ box, ncol = 1,
               labeller = labeller(box = labels_bio)) +
    ggtitle("Biomass") +
    theme_bw() + 
    theme(
      panel.grid = element_blank(),
      text = element_text(color = "black", size = 10),
      panel.border = element_rect(color = "black", linewidth = 0.5),
      legend.title = element_text(size = 10),
      aspect.ratio = 1,
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold", vjust = -5),
      strip.text = element_text(hjust = -0.03, vjust = 0, face = "bold", size = 10),
      panel.spacing.y = unit(0, "cm"),
      plot.margin = unit(c(-0.1,0.1,0,0), "lines"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9, color = "black"))
  
  exp_plot <- plot_a + plot_b

  fig_list <- c(3,4,7,10)
  ggsave(
    paste0(fig_dir, "Fig_", fig_list[g], ".tiff"),
    dpi = 600,
    exp_plot,
    width = 5,
    height = 6.5, 
    bg = "white")
}

### Diversity Plots ----
alpha_all <- do.call("rbind", lapply(1:30, function(x) {
  data.frame(
    lat = paste("Lat.",1:bx,sep=""),
    value = all_analyses[[x]]$div$alpha,
    all_analyses[[x]]$div$ids)
})) %>%
  mutate(t_amp = as.numeric(as.character(factor(lat, levels = unique(lat), labels = apply(
    temp_mat, 2, function(x){diff(range(x))}))))) %>%
  mutate(t_amp = ifelse(var == "T0", 0, t_amp))
alpha_all$lat<-factor(alpha_all$lat,levels = unique(alpha_all$lat))
alpha_all$factor.m<-factor(alpha_all$m,levels=unique(alpha_all$m),labels=mort_labs)
alpha_all$factor.m<-factor(alpha_all$factor.m,levels= levels(alpha_all$factor.m)[order(levels(alpha_all$factor.m))])


gamma_all <- do.call("rbind", lapply(1:30, function(x) {
  data.frame(
    lat = paste("Lat.",1:bx,sep=""),
    value = all_analyses[[x]]$div$gamma,
    all_analyses[[x]]$div$ids)
})) %>%
  mutate(t_amp = as.numeric(as.character(factor(lat, levels = unique(lat), labels = apply(
    temp_mat, 2, function(x){diff(range(x))}))))) %>%
  mutate(t_amp = ifelse(var == "T0", 0, t_amp))
gamma_all$lat<-factor(gamma_all$lat,levels = unique(gamma_all$lat))
gamma_all$factor.m<-factor(gamma_all$m,levels=unique(gamma_all$m),labels=mort_labs)
gamma_all$factor.m<-factor(gamma_all$factor.m,levels= levels(gamma_all$factor.m)[order(levels(gamma_all$factor.m))])

data_text <- data.frame(
  factor.m = mort,
  label = paste0(letters[1:length(mort)],".")
)
data_text$factor.m<-factor(data_text$factor.m,levels=unique(data_text$factor.m),labels=mort_labs)
data_text$factor.m<-factor(data_text$factor.m,levels=levels(data_text$factor.m)[order(mort)])


###### Fig. 5 - T0_D1 ------

plot_gdiv <- ggplot(alpha_all %>%
                      filter(var == "T0" & disp == "D1")) +
  geom_tile(aes(x = mag, y = lat, fill = value)) +
  scale_fill_gradientn(
    name = bquote(atop("5-year\naverage\nrichness",(bar(S)))),
    colours = brewer.pal(9,"Greens")[3:9],
    na.value = "grey",
    breaks = pretty(
      seq(1,
          max(as.numeric(as.character(alpha_all$value))),
          length.out = 5))) +
  scale_y_discrete(
    name = "Latitude",
    expand = c(0,0),
    labels = lats_plot,
    breaks = paste(
      "Lat",
      which(full_latitude%in%lats_plot),
      sep="."))+
  scale_x_discrete(
    name = bquote(K[H]~(m^2~s^-1)),
    expand = c(0,0),
    labels = c(expression(atop(10^0, bold("           Low dispersal"))),
               expression(10^1),
               expression(10^2),
               expression(atop(10^3, bold("High dispersal           ")))),
    breaks = paste("E", 1:length(k_h), sep="")) +
  theme_bw()+ 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold", vjust = -1),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black", vjust = 0)) +
  coord_cartesian(
    clip = "off",
    ylim = c(1,159))+
  geom_text(
    data = data_text,
    aes(x = 0.5, y = 164, label = label), 
    fontface = "bold",
    size = 4)+
  facet_grid(~factor.m, labeller = label_parsed)

ggsave(
  paste0(fig_dir, "Fig_5.tiff"),
  dpi=600,
  plot_gdiv,
  width=7.5,
  height=4.5,
  units = "in")

###### Fig. 8 - T1_D0 ------
plot_a <- ggplot(
  alpha_all %>% filter(var == "T1" & disp == "D0"),
  aes(
    x = t_amp,
    y = value,
    col = m,
    group = m))+
  geom_point() +
  geom_smooth(
    aes(group = m),
    method = "gam",
    formula = y~poly(x,3),
    se = F) +
  scale_y_continuous(name=bquote(atop(bold(bar(S)),5*'-'*year~average)), limits = c(1, NA))+
  scale_x_continuous(
    name = bquote(
      Temperature~amplitude~(degree*C)))+
  scale_color_manual(
    name = bquote(italic(gamma)),
    labels = c("exp1" = "0.05", "exp2" = "0.1", "exp3" = "0.2"),
    values = niche_colors)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    text = element_text(color="black",size=10),
    panel.border = element_rect(color="black"),
    axis.title = element_text(size=10),
    axis.text = element_text(size=9,color="black"),
    plot.margin = unit(c(1,0.5,0.5,0.5),"lines"),
    legend.position = "none")

plot_b <- ggplot(
  gamma_all %>% filter(var == "T1" & disp == "D0"),
  aes(
    x = t_amp,
    y = value,
    col = m,
    group = m))+
  geom_point()+
  geom_smooth(
    aes(group = m),
    method = "gam",
    formula = y ~ poly(x,3),
    se = F) +
  scale_y_continuous(name=bquote(atop(bold(S[T]),5*'-'*year~total)), expand = c(0.01,0.01))+
  scale_x_continuous(name=bquote(
    Temperature~amplitude~(degree*C)))+
  scale_color_manual(
    name = bquote(italic(gamma)),
    labels = c("exp1" = "0.05", "exp2" = "0.1", "exp3" = "0.2"),
    values = niche_colors)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    text = element_text(color="black",size=10),
    panel.border = element_rect(color="black"),
    axis.title = element_text(size=10),
    axis.text = element_text(size=9,color="black"),
    plot.margin = unit(c(1,0.5,0.5,0.5),"lines"),
    legend.position = c(0.2,0.88), 
    legend.background = element_rect("transparent"),
    legend.spacing.y = unit(0,"lines"),
    legend.title = element_text(size=9),
    legend.key.height = unit(0.1, "in"))

div_temp <- plot_b + plot_a & plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(face = 'bold'))

ggsave(
  paste0(fig_dir, "Fig_8.tiff"),
  dpi = 600,
  div_temp,
  width = 7,
  height = 3.5)


###### Fig. 11 - T1_D1 ------

alpha_all$div <- "atop(bold(bar(S)),5*'-'*year~average)"
gamma_all$div <- "atop(bold(S[T]),5*'-'*year~total)"
total_div <- rbind(alpha_all, gamma_all)
data_text <- data.frame(
  factor.m = rep(mort, 2),
  div = rep(unique(total_div$div), each=3),
  label = paste0(letters[1:6],".")
)
data_text$factor.m <- factor(
  data_text$factor.m,
  levels = unique(data_text$factor.m),
  labels = mort_labs)

combo_heat <- ggplot(total_div %>% filter(var == "T1" & disp == "D1"))+
  geom_tile(aes(
    x = mag,
    y = lat,
    fill = as.numeric(as.character(value))))+
  scale_fill_gradientn(
    name = "# Species",
    colours = brewer.pal(9,"Greens")[3:9],
    na.value = "grey",
    breaks = pretty(
      seq(1, max(
        as.numeric(
          as.character(alpha_all$value))),length.out = 5)
    ))+
  scale_y_discrete(
    name = bquote(Latitude),
    expand = c(0,0),
    labels = lats_plot,
    breaks = paste(
      "Lat",
      which(full_latitude%in%lats_plot),
      sep="."))+
  scale_x_discrete(
    name = bquote(K[H]~(m^2~s^-1)),
    expand = c(0,0),
    labels = c(expression(atop(10^0, bold("           Low dispersal"))),
               expression(10^1),
               expression(10^2),
               expression(atop(10^3, bold("High dispersal           ")))),
    breaks = paste("E", 1:length(k_h), sep="")) +
  theme_bw()+ 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold", vjust = -0.5),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black"),
    panel.spacing.y = unit(1,"lines")) +
  coord_cartesian(clip="off", ylim = c(1,159))+
  geom_text(data=data_text, aes(
    x = 0.5,
    y = 167,
    label = label),
    fontface = "bold",
    size = 4)+
  facet_grid(div~factor.m, labeller = label_parsed)

ggsave(
  paste0(fig_dir,"Fig_11.tiff"),
  dpi = 600,
  combo_heat,
  width = 7.5,
  height = 4.5)

### Niche Plots ----
for (i in 1:num_species) {
  converged_vals$t_rat[
    converged_vals$name == i] <- (converged_vals$tobs[converged_vals$name == i] -
                                    fund_topt[i]) / fund_topt[i]
  converged_vals$w_rat[
    converged_vals$name == i] <- (converged_vals$wobs[converged_vals$name == i] -
                                    fund_widths[i]) / fund_widths[i]
}
converged_vals$ogmag <- as.numeric(
  str_extract(converged_vals$mag, "\\d"))
converged_vals$mag <- as.numeric(
  str_extract(converged_vals$mag, "\\d"))

for (f in 1:length(mort)) {
  jitter <- seq(-0.2, 0.2, length.out = length(mort))[f]
  converged_vals$mag[
    converged_vals$m == paste0("exp",f)] <- converged_vals$mag[
      converged_vals$m == paste0("exp",f)] + jitter
}

converged_vals$factor.m <- factor(
  converged_vals$m,
  levels = unique(converged_vals$m)[order(mort)])

id_list <- do.call("rbind",(lapply(all_analyses,function(f){f$div$ids})))
exp_id <- which(apply(id_list, 1, function(r){all(c("T1","D0","E0") %in% r)}))
converged_vals$tdiff2 <- NA
for (l in 1:3){
  for (i in seq_along(unique(converged_vals$name))){
    box_pres <- which(unlist(lapply(all_analyses[[exp_id[l]]]$div$species, function(r){any(unique(converged_vals$name)[i] %in% r)})))
    converged_vals$tdiff2[converged_vals$name == unique(converged_vals$name)[i] & converged_vals$var == "T1" & converged_vals$disp == "D0" & converged_vals$m == paste0("exp",l)] <- diff(range(temp_mat[time_lastyr,box_pres], na.rm = TRUE))
  }
}
###### Fig. 6 - T0_D1 ------
fig_6_topt <- ggplot(
  data = converged_vals %>% filter(var == "T0" & disp == "D1"),
  aes(
    x = jitter(mag),
    y = t_rat,
    color = factor.m))+
  geom_point(alpha = 0.3)+
  geom_boxplot(aes(
    group = mag,
    fill = factor.m),
    color = "black",
    size = 0.2,
    outlier.color = NA,
    show.legend = F) +
  scale_y_continuous(name = bquote(delta[T^opt]))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_continuous(
    name = bquote(K[H]~(m^2~s^-1)),expand = c(0.005,0.005),
    labels = c(expression(atop(10^0, bold("           Low dispersal"))),
               expression(10^1),
               expression(10^2),
               expression(atop(10^3, bold("High dispersal           ")))),
    breaks = 1:4)+
  scale_color_manual(
    name = bquote(italic(gamma)),
    values = niche_colors,
    labels = as.numeric(mort)[order(mort)])+
  scale_fill_manual(
    name = bquote(italic(gamma)),
    values = niche_colors,
    labels = as.numeric(mort)[order(mort)])+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    text = element_text(color="black",size=10),
    panel.border = element_rect(color="black"),
    axis.title = element_text(size=10),
    axis.text = element_text(size=9,color="black"),
    plot.margin = unit(c(1,0.5,0.5,0.5),"lines"),
    legend.title = element_text(size = 12),
    legend.position = "none",
    axis.title.y = element_text(size = 13))



fig_6_w <- ggplot(
  data = converged_vals %>% filter(var == "T0" & disp == "D1"),
  aes(
    x = jitter(mag),
    y = w_rat,
    color = factor.m))+
  geom_point(alpha=0.3)+
  geom_boxplot(
    aes(group = mag, fill = factor.m),
    size = 0.2,
    color = "black",
    outlier.color = NA,
    show.legend = F)+
  geom_hline(
    yintercept = 0,
    linetype = "dashed")+
  scale_y_continuous(name = bquote(delta[W]))+
  scale_x_continuous(
    name = bquote(K[H]~(m^2~s^-1)),expand = c(0.005,0.005),
    labels = c(expression(atop(10^0, bold("           Low dispersal"))),
               expression(10^1),
               expression(10^2),
               expression(atop(10^3, bold("High dispersal           ")))),
    breaks = 1:4)+
  scale_color_manual(
    name = bquote(italic(gamma)),
    values = niche_colors,
    labels = as.numeric(mort)[order(mort)])+
  scale_fill_manual(
    name = bquote(italic(gamma)),
    values = niche_colors,
    labels = as.numeric(mort)[order(mort)])+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(color="black",size=10),
    panel.border = element_rect(color="black"),
    axis.title = element_text(size=10),
    axis.text = element_text(size=9,color="black"),
    plot.margin = unit(c(1,0.5,0.5,0.5),"lines"),
    legend.position = c(0.1,0.88), 
    axis.title.y = element_text(size=13),
    legend.background = element_rect("transparent"),
    legend.spacing.y = unit(0,"lines"),
    legend.spacing.x = unit(0, "lines"),
    legend.title = element_text(size = 11, hjust = 0.5),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.1, "in"))+
  guides(colour = guide_legend(
    override.aes = list(alpha = 1)))

fig_6_niche_disp <- ggarrange(
  fig_6_w, fig_6_topt, labels = c("a.","b."), ncol = 2, font.label = list(size = 10))

ggsave(
  paste0(
    fig_dir,
    "Fig_6.tiff"),
  dpi = 600,
  fig_6_niche_disp,
  width = 6.5,
  height = 3.5)

##### Fig. 9 - T1_D0 ----
fig_9_topt <- ggplot(
  data = converged_vals %>% filter(var == "T1" & disp == "D0"),
  aes(
    x = tdiff2,
    y = t_rat,
    col = factor.m))+
  geom_point()+
  scale_y_continuous(name = bquote(delta[T^opt]))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_continuous(
    name = bquote(
      Temperature~amplitude~(degree*C)))+
  scale_color_manual(
    name = bquote(italic(gamma)),
    values = niche_colors,
    labels = as.numeric(mort)[order(mort)])+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    text = element_text(color="black",size=10),
    panel.border = element_rect(color="black"),
    axis.title = element_text(size=10),
    axis.text = element_text(size=9,color="black"),
    plot.margin = unit(c(1,0.5,0.5,0.5),"lines"),
    legend.position = "none",
    axis.title.y = element_text(size=13))

fig_9_w <- ggplot(
  data = converged_vals %>% filter(var == "T1" & disp == "D0"),
  aes(
    x = tdiff2,
    y = w_rat,
    col = factor.m))+
  geom_point()+
  geom_smooth(
    method="gam",
    formula = y ~ poly(x,2),
    aes(group=factor.m),
    se=F)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_y_continuous(name = bquote(delta[W]))+
  scale_x_continuous(name = bquote(
    Temperature~amplitude~(degree*C)))+
  scale_color_manual(
    name = bquote(italic(gamma)),
    values = niche_colors,
    labels = as.numeric(mort)[order(mort)])+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    text = element_text(color="black",size=10),
    panel.border = element_rect(color="black"),
    axis.title = element_text(size=10),
    axis.text = element_text(size=9,color="black"),
    plot.margin = unit(c(1,0.5,0.5,0.5),"lines"),
    legend.position = c(0.2,0.88), 
    legend.background = element_rect("transparent"),
    legend.spacing.y = unit(0,"lines"),
    axis.title.y = element_text(size=13),
    legend.title = element_text(size=9, hjust = 0.5),
    legend.key.height = unit(0.1, "in"))

fig_9_niche_tamp <- ggarrange(
  fig_9_w, fig_9_topt, labels = c("a.","b."), ncol = 2, font.label = list(size = 10))

ggsave(
  paste0(fig_dir,"Fig_9.tiff"),
  dpi = 600,
  fig_9_niche_tamp,
  width = 6.5,
  height = 3)

##### Fig. 12 - T1_D1 -----

fig_12_w <- ggplot(
  data = converged_vals %>% filter(var == "T1" & disp == "D1"),
  aes(
    x = mag,
    y = w_rat,
    fill = tdiff,
    shape = m))+
  geom_point(
    color = "white",
    size = 4)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_fill_gradientn(
    name = bquote(Delta*T),
    colors = rev(brewer.pal(9, "Spectral")))+
  scale_y_continuous(name = bquote(delta[W]))+
  scale_x_continuous(
    name = bquote(K[H]~(m^2~s^-1)),expand = c(0.005,0.005),
    labels = c(expression(atop(10^0, bold("           Low dispersal"))),
               expression(10^1),
               expression(10^2),
               expression(atop(10^3, bold("High dispersal                 ")))),
    breaks = 1:4)+
  scale_shape_manual(
    name = bquote(italic(gamma)),
    labels = c("exp1" = "0.05", "exp2" = "0.1", "exp3" = "0.2"),
    values = c(21,22,24),
    guide = guide_legend(
      override.aes = list(color="black")))+
  scale_linetype_manual(
    name = bquote(italic(gamma)),
    labels = c("exp1" = "0.05", "exp2" = "0.1", "exp3" = "0.2"),
    values = c(1,3,6))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    text = element_text(color="black",size=11),
    panel.border = element_rect(color="black"),
    axis.title = element_text(size=11),
    axis.text = element_text(size=10,color="black"),
    plot.margin = unit(c(1,1,0.5,0.5),"lines"),
    legend.title = element_text(size = 13, hjust = 0.5),
    legend.spacing.x = unit(0.2, "lines"),
    axis.title.y = element_text(size=13))


fig_12_topt <- ggplot(
  data = converged_vals %>% filter(var == "T1" & disp == "D1"),
  aes(
    x = mag,
    y = t_rat,
    fill = tdiff,
    shape = m))+
  geom_point(color = "white", size = 4)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_fill_gradientn(
    name = bquote(Delta*T),
    colors = rev(brewer.pal(9, "Spectral")))+
  scale_y_continuous(name = bquote(delta[T^opt]))+
  scale_x_continuous(
    name = bquote(K[H]~(m^2~s^-1)),expand = c(0.005,0.005),
    labels = c(expression(atop(10^0, bold("           Low dispersal"))),
               expression(10^1),
               expression(10^2),
               expression(atop(10^3, bold("High dispersal           ")))),
    breaks = 1:4)+
  scale_shape_manual(
    name = bquote(italic(gamma)),
    labels = c("exp1" = "0.05", "exp2" = "0.1", "exp3" = "0.2"),
    values = c(21,22,24),
    guide = guide_legend(
      override.aes = list(color = "black")))+
  scale_linetype_manual(
    name = bquote(italic(gamma)),
    labels = c("exp1" = "0.05", "exp2" = "0.1", "exp3" = "0.2"),
    values = c(1,3,6))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    text = element_text(color="black",size=11),
    panel.border = element_rect(color="black"),
    axis.title = element_text(size=11),
    axis.text = element_text(size=10,color="black"),
    plot.margin = unit(c(1,0.5,1,0.5),"lines"),
    legend.title = element_text(size = 13, hjust = 0.5),
    legend.text.align = 0,
    axis.title.y = element_text(size=13),
    legend.position = "none")


fig_12_niche_tamp_disp <- fig_12_w + plot_spacer() +fig_12_topt + plot_layout(guides = 'collect', widths = c(1,0.1,1)) & plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(face = 'bold'), plot.margin = margin(0.1,0.1,0.1,0.1, "lines"))

ggsave(
  paste0(fig_dir,"Fig_12.tiff"),
  dpi = 320,
  fig_12_niche_tamp_disp,
  width = 7,
  height = 3.5)


### Fig. 13 - Summary figure ----

summ_a <- subset(converged_vals, m == "exp1") %>%
  mutate(fac = factor(
    paste0(var, disp, m),
    levels = c("T0D0exp1", "T0D1exp1", "T1D0exp1", "T1D1exp1"),
    labels = c("E1", "E2", "E3", "E4"))) %>%
  ggplot() + 
  geom_point(
    aes(x = fac, y = w_rat, color = factor(ogmag),
        alpha = factor(ogmag)),
    position = position_jitterdodge()) +
  geom_boxplot(
    aes(x = fac, y = w_rat, fill = factor(ogmag)),
    outlier.colour = NA,
    color = "black") +
  labs(x = "Experiment", y = bquote(delta[W])) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(color = "black", size = 13),
        panel.grid = element_blank()) +
  scale_fill_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    labels = c(
      "0" = "None",
      "1" = bquote(10^0),
      "2" = bquote(10^1),
      "3" = bquote(10^2),
      "4" = bquote(10^3)),
    values = c("white", brewer.pal(5, "Blues")[2:5])
    )+
  scale_color_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    labels = c(
      "0" = "None",
      "1" = bquote(10^0),
      "2" = bquote(10^1),
      "3" = bquote(10^2),
      "4" = bquote(10^3)), 
    values = c("black", brewer.pal(5, "Blues")[2:5]),
    guide = "none") +
  scale_alpha_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    values = c(
      "0" = 0.3,
      "1" = 0.8,
      "2" = 0.8,
      "3" = 0.8,
      "4" = 0.8), 
    guide = "none") +
  theme(legend.position = c(0.15,0.75),
        legend.text.align = 0.1,
        legend.title = element_text(size = 10),
        legend.background = element_blank())

summ_b <- subset(converged_vals, m == "exp1") %>%
  mutate(fac = factor(
    paste0(var, disp, m),
    levels = c("T0D0exp1", "T0D1exp1", "T1D0exp1", "T1D1exp1"),
    labels = c("E1", "E2", "E3", "E4"))) %>%
  ggplot() + 
  geom_point(
    aes(x = fac, y = t_rat, color = factor(ogmag), alpha = factor(ogmag)),
    position = position_jitterdodge()) +
  geom_boxplot(
    aes(x = fac, y = t_rat, fill = factor(ogmag)),
    outlier.colour = NA,
    color = "black") +
  labs(x = "Experiment", y = bquote(delta[T^opt])) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(color = "black", size = 13),
        panel.grid = element_blank()) +
  scale_fill_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    labels = c(
      "0" = "None",
      "1" = bquote(10^0),
      "2" = bquote(10^1),
      "3" = bquote(10^2),
      "4" = bquote(10^3)),
    values = c("white", brewer.pal(5, "Blues")[2:5])
  )+
  scale_color_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    labels = c(
      "0" = "None",
      "1" = bquote(10^0),
      "2" = bquote(10^1),
      "3" = bquote(10^2),
      "4" = bquote(10^3)), 
    values = c("black", brewer.pal(5, "Blues")[2:5]),
    guide = "none") +
  scale_alpha_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    values = c(
      "0" = 0.3,
      "1" = 0.8,
      "2" = 0.8,
      "3" = 0.8,
      "4" = 0.8), 
    guide = "none") +
  theme(legend.position = "none")


summ_c <- subset(gamma_all, m == "exp1") %>%
  mutate(fac = factor(
    paste0(var, disp, m),
    levels = c("T0D0exp1", "T0D1exp1", "T1D0exp1", "T1D1exp1"),
    labels = c("E1", "E2", "E3", "E4"))) %>%
  ggplot() + 
  geom_point(
    aes(x = fac, y = value, color = factor(mag), alpha = factor(mag)),
    position = position_jitterdodge()) +
  geom_boxplot(
    aes(x = fac, y = value, fill = factor(mag)),
    outlier.colour = NA,
    color = "black") +
  labs(x = "Experiment",
       y = bquote(atop(S[T], Total~diversity~within~last~5~years))) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(color = "black", size = 13),
        panel.grid = element_blank()) +
  scale_fill_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    labels = c(
      "E0" = "None",
      "E1" = bquote(10^0),
      "E2" = bquote(10^1),
      "E3" = bquote(10^2),
      "E4" = bquote(10^3)),
    values = c("white", brewer.pal(5, "Blues")[2:5])
  )+
  scale_color_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    labels = c(
      "E0" = "None",
      "E1" = bquote(10^0),
      "E2" = bquote(10^1),
      "E3" = bquote(10^2),
      "E4" = bquote(10^3)), 
    values = c("black", brewer.pal(5, "Blues")[2:5]),
    guide = "none") +
  scale_alpha_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    values = c(
      "E0" = 0.3,
      "E1" = 0.8,
      "E2" = 0.8,
      "E3" = 0.8,
      "E4" = 0.8), 
    guide = "none") +
  theme(legend.position = "none")

summ_d <- subset(alpha_all, m == "exp1") %>%
  mutate(fac = factor(
    paste0(var, disp, m),
    levels = c("T0D0exp1", "T0D1exp1", "T1D0exp1", "T1D1exp1"),
    labels = c("E1", "E2", "E3", "E4"))) %>%
  ggplot() + 
  geom_point(
    aes(x = fac, y = value, color = factor(mag), alpha = factor(mag)),
    position = position_jitterdodge()) +
  geom_boxplot(
    aes(x = fac, y = value, fill = factor(mag)),
    outlier.colour = NA,
    color = "black") +
  labs(x = "Experiment",
       y = bquote(atop(bar(S), Avg~diversity~within~last~5~years))) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(color = "black", size = 13),
        panel.grid = element_blank()) +
  scale_fill_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    labels = c(
      "E0" = "None",
      "E1" = bquote(10^0),
      "E2" = bquote(10^1),
      "E3" = bquote(10^2),
      "E4" = bquote(10^3)),
    values = c("white", brewer.pal(5, "Blues")[2:5])
  )+
  scale_color_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    labels = c(
      "E0" = "None",
      "E1" = bquote(10^0),
      "E2" = bquote(10^1),
      "E3" = bquote(10^2),
      "E4" = bquote(10^3)), 
    values = c("black", brewer.pal(5, "Blues")[2:5]),
    guide = "none") +
  scale_alpha_manual(
    name = bquote(K[H]~(m^2~s^-1)),
    values = c(
      "E0" = 0.3,
      "E1" = 0.8,
      "E2" = 0.8,
      "E3" = 0.8,
      "E4" = 0.8), 
    guide = "none") +
  theme(legend.position = "none")

summary_box <- plot_grid(NULL, NULL, summ_a,summ_b,summ_c,summ_d, nrow = 3, labels = c(NA, NA, "a.", "b.", "c.","d."), label_y = 1.05, rel_heights = c(0.1,1,1))

ggsave(
  paste0(fig_dir, "Fig_13.tiff"),
  dpi = 600,
  summary_box,
  width = 8,
  height = 8,
  bg = "white")

### Fig. 14 - Source/sink dynamics -----
div_compare <- gamma_all %>%
  filter(m == "exp1" & mag %in% c("E0", "E3") & !(disp == "D0" & var == "T0")) %>%
  ggplot() +
  geom_point(aes(y = as.numeric(value), x = lat, color = paste(var, disp, mag, m, sep = "_")), show.legend = FALSE) +
  geom_smooth(
    aes(y = as.numeric(value), x = as.numeric(lat), color = paste(var, disp, mag, m, sep = "_")),
    method = "gam",
    formula = y ~ s(x, k = 6),
    fullrange = TRUE,
    se = F) +
  scale_x_discrete(
    name = "Latitude",
    expand = c(0,0),
    labels = lats_plot,
    breaks = paste(
      "Lat",
      which(full_latitude%in%lats_plot),
      sep=".")) +
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold", vjust = -1),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black")) +
  coord_flip(expand = 0,
             clip = "on", ylim = c(0.9,NA)) +
  labs(y = bquote(S[T])) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_color_manual(
    name = "Experiment",
    values = c(
      "T1_D1_E3_exp1" = "#F27059",
      "T0_D1_E3_exp1" = "#73BA9B",
      "T1_D0_E0_exp1" = "#FFD275"),
    labels = c(
      "T1_D1_E3_exp1" = "Combined",
      "T0_D1_E3_exp1" = "Dispersal",
      "T1_D0_E0_exp1" = "Temp. Variability"
    ))

niche_width_plot <- ggplot() +
  geom_smooth(
    data = converged_df_sp2 %>%
      mutate(temp = round(temp,2)) %>%
      group_by(temp, id) %>%
      summarise(bio = max(bio)) %>%
      mutate(id = factor(id, levels = c("T0_D1_E3_exp1", "T1_D0_E0_exp1",  "T1_D1_E3_exp1"))) %>% filter(id %in% c("T0_D1_E3_exp1", "T1_D1_E3_exp1")),
    aes(x = temp, y = bio, group = id, color = id),
    method = "gam", formula = y ~ poly(x,2), se = F) +
  geom_smooth(
    data = converged_df_sp2 %>%
      mutate(temp = round(temp,2)) %>%
      group_by(temp, id) %>%
      summarise(bio = max(bio)) %>%
      mutate(id = factor(id, levels = c("T0_D1_E3_exp1", "T1_D0_E0_exp1",  "T1_D1_E3_exp1"))) %>% filter(id == "T1_D0_E0_exp1"),
    aes(x = temp, y = bio, group = id, color = id),
    method = "gam", formula = y ~ poly(x,2), se = F, fullrange = T) +
  coord_cartesian(ylim = c(0, NA), expand = 0) +
  labs(x = bquote(Temperature~(degree*C)), y = bquote(Biomass~(mmolP~m^-3))) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.18, 0.85),
    legend.direction = "vertical",
    legend.background = element_blank()) +
  scale_color_manual(
    name = "",
    values = c(
      "T0_D1_E3_exp1" = "#73BA9B",
      "T1_D0_E0_exp1" = "#FFD275",
      "T1_D1_E3_exp1" = "#F27059"), #"#893168",
    labels = c(
      "T0_D1_E3_exp1" = "Experiment 2",
      "T1_D0_E0_exp1" = "Experiment 3",
      "T1_D1_E3_exp1" = "Experiment 4"
    ))

## Processes --
box_choose <- p_sink_df %>%
  group_by(box) %>%
  summarise(max = max(t_mu), time_max = time_vec[which.max(t_mu)]) %>%
  filter(box >= 80 & max > 0) %>% 
  mutate(max = round(max, 4)) %>% 
  summarise(
    box_happy = box[max %in% max(max)],
    time_happy = time_max[max %in% max(max)]) %>%
  slice(c(1,n())) %>%
  mutate(date = as.Date(365*(time_happy-4), origin = "2016-01-01"))

happy_sad <- p_sink_df %>%
  filter(time_vec %in% box_choose$time_happy) %>%
  mutate(
    trans = (trans*1200),
    bio = bio*1.5,
    time_lab = factor(
      time_vec,
      levels = box_choose$time_happy,
      labels = c(
        bquote(Maximum~growth~at~.(full_latitude[box_choose$box_happy[1]])*degree~N),
        bquote(Maximum~growth~at~.(full_latitude[box_choose$box_happy[2]])*degree~N)))) %>%
  select(-c(net_growth, growth, decay, nutr, temp)) %>%
  pivot_longer(!c(box,time_vec, time_lab), names_to = "var", values_to = "vals") %>%
  filter(box >= 80) %>%
  ggplot() +
  geom_line(aes(x = box, y = vals, color = var, linetype = ifelse(var == "t_mort", "solid", "dashed")), size = 0.8) +
  facet_wrap(~time_lab, labeller = label_parsed) +
  scale_y_continuous(name = bquote(days^-1),
                     limits = c(-3, 3),
                     sec.axis = sec_axis(trans=~./1.5, bquote(mmol*P~m^-3))) +
  scale_x_continuous(name = "Latitude", expand = c(0,0),
                     labels = lats_plot,
                     breaks = which(full_latitude%in%lats_plot)) +
  theme_bw() +
  scale_color_manual(
    name = "",
    values = c(
      "t_mu" = "black",
      "t_mort" = "black",
      "trans" = "#CC2936",
      "bio" = "#009DDC"
    ),
    labels = c(
      "t_mu" = "Growth",
      "trans" = "Net Transport",
      "t_mort" = "Mortality",
      "bio" = "Biomass")) +
  scale_linetype_discrete(guide = "none") +
  coord_flip() +
  theme(
    axis.text.x.bottom = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid = element_blank(),
    legend.position = c(0.12,0.84),
    legend.box = "veritcal",
    legend.background = element_blank(),
    legend.spacing.y = unit(0, "cm"),
    legend.key.width = unit(1,"cm"),
    legend.title = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.direction = "vertical",
    axis.line.x.top = element_line(color = "#009DDC"),
    axis.text.x.top = element_text(color =  "#009DDC"),
    axis.title.x.top = element_text(color =  "#009DDC")) +
  guides(color=guide_legend(nrow=4, byrow=TRUE, override.aes = list(linetype = c(1, 2, 1, 1), size = 0.5)))


final_compare_fig2 <- ((div_compare + niche_width_plot)/happy_sad & theme(axis.title = element_text(size = 10, color = "black"), axis.text = element_text(size = 9, color = "black")))

ggsave(paste0(getwd(), "/Figures/final_compare2-b.svg"), final_compare_fig2, width = 7, height = 8, units = "in")

### Supplementary Figures ----
# scp ans132@kaweah.ucsd.edu:documents/metacomm_model/rds.output/T1_D1_E2_exp1-publication.Rdata ~/Documents/PhD_Work/Research/Niche_Model/Data/
df_id <- c("T0_D0_E0_exp1", "T0_D1_E3_exp1", "T1_D0_E0_exp1", "T1_D1_E3_exp1")
id_list <- do.call("rbind",(lapply(all_analyses,function(f){f$div$ids})))

for (e in 1:4){
  # open nc file
  load(paste0(files_home,"rds.output/", df_id[e], "-publication.Rdata"))
  # identify which list element matches the variables
  exp_id <- which(apply(id_list, 1, function(r){all(str_split(df_id[e], "_")[[1]] %in% r)}))
  
  
  # Combined the time, temperature, and biomass for each surviving species
  # in each box into one dataframe
  plot_df <- do.call("rbind",lapply(seq_along(boxes_to_plot), function(r){
    df <- data.frame(
      time_vec = rep(1:iter, num_species),
      bio = matrix(P[,,boxes_to_plot[r]]),
      sp = rep(1:num_species, each = iter),
      box = boxes_to_plot[[r]])
    return(df)
  }))
  
  # Pull out the fundamental optimal temperatures for each unique species
  sp_use_plot <- fund_topt
  names(sp_use_plot) <- 1:num_species
  
  plot_df$box <- factor(plot_df$box, levels = rev(boxes_to_plot))
  # Set up tag labels (a., b., c., etc.)
  labels_bio <- paste0(letters[c(1:4)], ".")
  names(labels_bio) <- rev(boxes_to_plot)
  
  plot_df$sp <- factor(plot_df$sp, levels = unique(plot_df$sp)[order(unique(plot_df$sp))])
  sp_colors <- colors_P(length(unique(plot_df$sp)))
  sp_use_plot <- fund_topt[as.numeric(levels(plot_df$sp))]
  names(sp_use_plot) <- levels(plot_df$sp)
  names(sp_colors) <- levels(plot_df$sp)
  
  intro_label_temp <- data.frame(
    box = boxes_to_plot,
    label = paste0(full_latitude[boxes_to_plot], "°N")
  )
  
  plot_supp <- ggplot(plot_df) +
    geom_line(aes(
      x = time_vec*t_step/365,
      y = log10(bio),
      col = factor(sp),
      group = factor(sp))) +
    scale_color_manual(
      name = bquote(italic(T[F]^opt) ~ (degree * C)),
      labels = signif(sp_use_plot,3),
      values = sp_colors,
      guide = guide_legend(reverse = TRUE, ncol = 1, keyheight=0.3)) +
     scale_y_continuous(
      name = bquote(Biomass~(log[10]~mmol~P~m^-3))) +
    scale_x_continuous(
      name = bquote(Year),
      expand = c(0,0),
      breaks = c(1, 10, 20, 30, 40, 50)) +
    coord_cartesian(ylim = c(-6, 1.5), expand = F) +
    facet_wrap(. ~ box, ncol = 1,
               labeller = labeller(box = labels_bio)) +
    ggtitle("Biomass") +
    theme_bw() + 
    theme(
      panel.grid = element_blank(),
      text = element_text(color = "black", size = 10),
      panel.border = element_rect(color = "black", linewidth = 0.5),
      legend.title = element_text(size = 10),
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold", vjust = -5),
      strip.text = element_text(hjust = 0, vjust = 0, face = "bold", size = 10),
      panel.spacing.y = unit(0, "cm"),
      plot.margin = unit(c(-0.1,0.1,0,0), "lines"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9, color = "black"))
  
  ggsave(
    paste0(files_home, "supplement_", e, "-full_times.tiff"),
    dpi = 600,
    plot_supp,
    width = 6.5,
    height = 6.5, 
    bg = "white")
  
}

# scp ans132@kaweah.ucsd.edu:documents/metacomm_model/supplement_3-full_times.tiff ~/Documents/PhD_Work/Research/Niche_Model/Figures/Reviews/

# 
# Pau
# 
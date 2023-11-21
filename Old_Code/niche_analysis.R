#### Streamlined Figure Code for ANY Model Runs ----
#### Packages ----
packages <- c(
    "ncdf4", "tidyr", "class", "fields", "R.matlab",
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

### Functions ---
source("documents/metacomm_model/tniche_analysis-update.R")
source("documents/metacomm_model/diversity_function-update.R")


#### Folders ----
files_home <- c("documents/metacomm_model/ncdf/")
#files_home <- c("documents/rds.output/final/rerun/")

#### Metadata ----
pattern_search <- c("exp1", "exp2", "exp3")
mort <- c(1 / 20, 1 / 10, 1 / 5)
names(mort) <- paste("m_", seq_along(mort), sep = "")
mort_labs <- do.call("expression", lapply(
  mort, function(c) {
    substitute(italic(gamma) == X, list(X = c))}))
labs_kh <- do.call("expression", lapply(0:3, function(r) {
  substitute(10^X, list(X = r))}))
sp_sep <- 1
z_vals <- seq(-4, 40, by = sp_sep)
num_species <- length(z_vals) #species
t_step <- 3 / 24 #day
k <- 1 * (1 / 1E-3) * (1 / 1E3) #mmolPm^-3 half-saturation constant
d <- (1E-5) * 60 * 60 * 24 #vertical flux of nutrients from the deep; day^-1
R_0 <- (0.8) * (1 / 1E3) * (1 / 1E-3) #deep nutrients; chemostatic mmolP m^-3
years <- 50 #How many years the model runs for; should up to 100 on final runs
n_day <- 365 * years #Number of days for model run purposes
iter <- n_day / t_step #Model iterations
spacing <- 1
full_latitude <- seq(-79, 79, by = spacing) ### Discretation of full_latitudes
latitude <- seq(-79, 79 + spacing, by = spacing) #add one to round it out
bx <- length(full_latitude)


boxes_to_plot <- which(full_latitude %in% c(10, 35, 60, 75))
k_h <- c(c(1E0, 1E1, 1E2, 1E3)) * (60 * 60 * 24)

#### subtitle
nc_path <- files_home
subtitle <- "exp" #"50yrs_again"
file_list <- list.files(nc_path, pattern = subtitle)
file <- file_list[1]
time_vec <- seq(1, n_day, len = iter) #time_vec vector
oneyr_iter <- 365 / t_step
oneyear <- oneyr_iter
lastyears <- (iter - oneyr_iter * 5):iter
lastyear <- (iter - oneyr_iter):iter
yr_45 <- 365*45/t_step
time_lastyr <- (length(yr_45:iter) - oneyr_iter):length(yr_45:iter)
five_years <- 1:(oneyr_iter*5)

#### Fundamental Niche ----
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
fund_widths <- unlist(lapply(fund_vals, function(x) x$vals$width), use.names = F)
fund_topt <- unlist(lapply(fund_vals, function(x) x$vals$topt), use.names = F)
fund_curves <- do.call("rbind", lapply(fund_vals, function(x) x$df))

### Run analysis ---
ensemble_files <- list.files(files_home)[str_detect(
  list.files(files_home), "ensemble")]
ensemble_analysis <- lapply(ensemble_files, tniche_analysis)


save(control_analysis, file = paste0(
  files_home, "control_analysis", subtitle, ".Rdata"))



#### E1/Control: No temperature Analysis ----
  control_files <- file_list[str_detect(file_list, paste0("T0_D0_"))]
  control_analysis <- lapply(control_files, tniche_analysis)
  save(control_analysis, file = paste0(
    files_home, "control_analysis", subtitle, ".Rdata"))

#### E2: Spatial Mass Effects Analysis ----
  spatial_files <- list.files(files_home)[
    str_detect(list.files(files_home), paste0("T0_D1_"))]
  spatial_analysis <- lapply(spatial_files, tniche_analysis)
  save(spatial_analysis, file = paste0(
    files_home, "spatial_analysis-", subtitle, ".Rdata"))
#### E3: Temporal Mass Effect Analysis ----
  temporal_files <- list.files(files_home)[str_detect(list.files(files_home), paste0("T1_D0_"))]
  temporal_analysis <- lapply(temporal_files, tniche_analysis)
  save(temporal_analysis, file = paste0(
    files_home, "temporal_analysis", subtitle, ".Rdata"))
### E4: Combined Effects Analysis ----
  combined_files <- file_list[grepl("T1_D1_", file_list) & !grepl("year30", file_list)]
  combined_analysis <- lapply(combined_files, tniche_analysis)
  save(combined_analysis, file = paste0(
    files_home, "combined_analysis-", subtitle, ".Rdata"))
  
  ensemble_files <- paste0(files_home, "T0_D0_E0_exp1-review_ensemble.nc")

  
### All together
  files_exp3 <-list_files(files_home)[grep("*D1_E[1-4]_exp3", list_files(files_home))][-c(1,3)]
  files_fixed <- list_files(files_home)[grep("D0_E0_exp[1-3]-publication-2", list_files(files_home))]
  files_rest <- list_files(files_home)[grep("*D1_E[1-4]_exp[1-2]-publication-2", list_files(files_home))]
  list_all2 <- c(files_fixed,files_rest,files_exp3)
  all_analyses <- lapply(list_all2, tniche_analysis)
  save(all_analyses, file = "./documents/metacomm_model/csv_output/all_diversity_list-publication-2.Rdata")
  

  file <- c("T0_D0_E0_exp2-publication-2.nc")
  testing_ext <- lapply(file, tniche_analysis)
  save(testing_ext, file = "./documents/metacomm_model/csv_output/all_diversity_list-publication-mini.Rdata")
    ## For comparison figure
    files_compare <- c(1, 10, 16, 25)
  
#### Load ----
lapply(list(
  "combined_analysis",
  "temporal_analysis",
  "spatial_analysis",
  "control_analysis"), function(u) {
    load(
      paste0(files_home, list_files(
        files_home, pattern = paste0(u, "-", subtitle, ".Rdata"))),
        .GlobalEnv)})

#### E1: Control Figures ----
only_t <- str_detect(file_list, "T0_D0_E0_exp1")
load(paste(files_home, file_list[only_t], sep = ""))

p_timeavg <- array(apply(P[1:iter,,], 3, function(x) {
  apply(x, 2, function(y) {
    as.vector(by(y, ceiling(1:iter / oneyr_iter), function(z) {
      mean(z, na.rm = T)}))})}), dim=c(years, n, bx))
p_timediff<-array(apply(signif(p_timeavg,2),3,diff),dim=c(years-1,n,bx))
div.no<-diversity_function(x=p_timediff[(years-5):(years-1),,],y=p_timeavg,z=P[lastyear,,])



l<-boxes_to_plot
p1<-do.call("rbind",lapply(1:length(l),function(c) {
  if(length(unique(unlist(div.no[[1]][[l[c]]])))==1) {
    df<-data.frame(box=l[c],time_vec=as.vector(time_vec[lastyear]/365),temp=as.vector(t[lastyear,l[c]]),P[lastyear,unlist(div.no[[1]][[l[c]]]),l[c]])
    colnames(df)[4]<-paste0("Species.",unlist(div.no[[1]][[l[c]]])) 
  }else{
    df<-data.frame(box=l[c],time_vec=as.vector(time_vec[lastyear]/365),temp=as.vector(TEMP[lastyear,l[c]]),P[lastyear,unlist(div.no[[1]][[l[c]]]),l[c]])
  }
  df%>%pivot_longer(cols = !c("box","time_vec","temp"))
}))
plot.df<-p1
plot.df$name<-factor(plot.df$name,levels=paste("Species",1:n,sep="."))
labels.temp<-paste(letters[seq(1,length(l)*2,by=2)],".",sep="")
names(labels.temp)<-l
labels.ab<-paste(letters[seq(2,(length(l)*2),by=2)],".",sep="")
names(labels.ab)<-l

poc.par<-theme(panel.grid = element_blank(),
                text = element_text(color="black",size=10),
                panel.border = element_rect(color="black"),
                legend.title = element_text(size=10),
                aspect.ratio=1,
                strip.background = element_blank(),
                strip.text = element_text(hjust=0,face="bold",size=10),
                panel.spacing.y = unit(0,"cm"),
                plot.margin = rep(unit(0,"null"),4),
                axis.title = element_text(size=10),
                axis.text = element_text(size=9,color="black"))

plot.a<-ggplot(plot.df)+
  geom_line(aes(x=time_vec,y=temp))+
  scale_y_continuous(name=bquote(Temperature~(degree*C)))+
  scale_x_continuous(name=bquote(time_vec~(year)),expand=c(0,0),breaks=quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),labels=signif(quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),3))+
  facet_wrap(~box,ncol = 1,labeller = labeller(box=labels.temp))+
  theme_bw()+poc.par
  
sp.names<-unique(plot.df$name)[order(unique(plot.df$name))]
sp.use.plot<-gsub("Species.","",sp.names)
sp.use.plot.2<-z_vals[as.numeric(sp.use.plot)][order(z_vals[as.numeric(sp.use.plot)])]
names(sp.use.plot.2)<-sp.names[order(as.numeric(sp.use.plot))]
sp.names.none<-sp.names

plot.b<-ggplot(plot.df)+
  geom_line(aes(x=time_vec,y=log10(value),col=name,group=name))+
  scale_color_manual(name=bquote(T[opt]~(degree*C)),
                      labels=sp.use.plot.2,values=colors_P(length(sp.use.plot.2)))+
  scale_y_continuous(name=bquote(log[10]~Biomass~(mmol~P~m^-3)),expand=c(0,0), limits=c(-4,1.5))+
  scale_x_continuous(name=bquote(time_vec~(year)),expand=c(0,0),breaks=quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),labels=signif(quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),3))+
  facet_wrap(~box,ncol = 1,labeller = labeller(box=labels.ab))+
  theme_bw()+poc.par

exp0=list(a=plot.a,b=plot.b)
save(exp0,file=paste0(proof_files,"exp0_plot-",subtitle,".Rdata"))


load(paste0(files_home,"exp0_plot-",subtitle,".Rdata"))
proof.none<-ggarrange(exp0$a,exp0$b,ncol=2,common.legend = T,legend="right")
ggsave(paste(figure_home,"proof_none-",subtitle,".svg",sep=""),dpi=600,proof.none,width=4,height=6)

### Summary Values ----
alpha.no<-do.call("rbind",lapply(1:(length(mort)),function(x) {data.frame(lat=paste("Lat.",1:bx,sep=""),value=control_analysis[[x]]$div$alpha,control_analysis[[x]]$div$ids)}))


#### E2: Spatial Mass Effects Figures ----
  ### Proof of Concept Plot -----
  only_t<-str_detect(file_list, "T0_D1_E1_exp1")
  load(paste(files_home,file_list[only_t],sep=""))

  div.spa<-spatial_analysis[[1]][[2]]$species
  
  l<-boxes_to_plot
  p1<-do.call("rbind",lapply(1:length(l),function(c) {
    data.frame(box=l[c],time_vec=as.vector(time_vec[lastyear]/365),temp=as.vector(TEMP[lastyear,l[c]]),P[lastyear,unlist(div.spa[[l[c]]]),l[c]])%>%
      pivot_longer(cols = !c("box","time_vec","temp"))
  }))
  #p1%>%group_by(box)%>%summarise(unique(name))
  plot.df<-p1
  plot.df$name<-factor(plot.df$name,levels=paste("Species",1:n,sep="."))
  plot.df$lt<-"np"
  plot.df$lt[plot.df$name %in% sp.names.none & plot.df$value>1E-3]="p"
  labels.temp<-paste(letters[seq(1,length(l)*2,by=2)],".",sep="")
  names(labels.temp)<-l
  labels.ab<-paste(letters[seq(2,(length(l)*2),by=2)],".",sep="")
  names(labels.ab)<-l
  
  
  plot.a<-ggplot(plot.df)+
    geom_line(aes(x=time_vec,y=temp))+
    scale_y_continuous(name=bquote(Temperature~(degree*C)))+
    scale_x_continuous(name=bquote(time_vec~(year)),expand=c(0,0),breaks=quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),labels=signif(quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),3))+
    facet_wrap(~box,ncol = 1,labeller = labeller(box=labels.temp))+
    theme_bw()+poc.par
  sp.names<-unique(plot.df$name)[order(unique(plot.df$name))]
  sp.use.plot<-gsub("Species.","",sp.names)
  sp.use.plot.2<-z_vals[as.numeric(sp.use.plot)][order(z_vals[as.numeric(sp.use.plot)])]
  names(sp.use.plot.2)<-sp.names[order(as.numeric(sp.use.plot))]
  plot.b<-ggplot(plot.df)+
    geom_line(aes(x=time_vec,y=log10(value),col=name,group=name,linetype=lt))+
    scale_color_manual(name=bquote(T[opt]~(degree*C)),
                       labels=sp.use.plot.2,values=colors_P(length(sp.use.plot.2)))+
    scale_linetype_manual(values=c("np"="solid","p"="dotdash"),guide=F)+
    scale_y_continuous(name=bquote(log[10]~Biomass~(mmol~P~m^-3)),expand=c(0,0), limits=c(-4,1.5))+
    scale_x_continuous(name=bquote(time_vec~(year)),expand=c(0,0),breaks=quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),labels=signif(quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),3))+
    facet_wrap(~box,ncol = 1,labeller = labeller(box=labels.ab))+
    theme_bw()+poc.par
  
  
  exp_spatial=list(a=plot.a,b=plot.b)
  save(exp_spatial,file=paste0(proof_files,"exp-spatial_plot-",subtitle,".Rdata"))
  
  
  load(paste0(files_home,"exp-spatial_plot-",subtitle,".Rdata"))
  proof.disp<-ggarrange(exp_spatial$a,exp_spatial$b,ncol=2,common.legend = T,legend="right")
  ggsave(paste(figure_home,"proof_disp.svg",sep=""),dpi=600,proof.disp,width=4,height=6)
  
  ### Diversity Plot ----
  alpha.all<-do.call("rbind",lapply(1:(length(k_h)*length(mort)),function(x) {data.frame(lat=paste("Lat.",1:bx,sep=""),value=spatial_analysis[[x]]$div$alpha,spatial_analysis[[x]]$div$ids)}))
  
  alpha.all$lat<-factor(alpha.all$lat,levels = unique(alpha.all$lat))
  alpha.all$factor.m<-factor(alpha.all$m,levels=unique(alpha.all$m),labels=mort_labs)
  
  alpha.all$factor.m<-factor(alpha.all$factor.m,levels= levels(alpha.all$factor.m)[order(levels(alpha.all$factor.m))])
  
  data_text <- data.frame(
    factor.m = mort,
    label = paste0(letters[1:length(mort)],".")
  )
  data_text$factor.m<-factor(data_text$factor.m,levels=unique(data_text$factor.m),labels=mort_labs)
  data_text$factor.m<-factor(data_text$factor.m,levels=levels(data_text$factor.m)[order(mort)])
  
  ###Stats
  spatial.glm<-gam(data=alpha.all,formula=value~mag*m)
  summary(spatial.glm)
 
plot.gdiv<-ggplot(alpha.all)+
    geom_tile(aes(x=mag,y=lat,fill=value))+
    scale_fill_gradientn(name=bquote(bar(S)~diversity),colours=brewer.pal(9,"Greens")[2:9],na.value = "grey",#limits=c(0.01,NA),
                         breaks=pretty(seq(1,max(as.numeric(as.character(alpha.all$value))),length_out=5)))+
    scale_y_discrete(name=bquote(Latitude),expand=c(0,0),
                     labels=lats_plot,
                     breaks=paste("Lat",which(full_latitude%in%lats_plot),sep="."))+
    scale_x_discrete(name=bquote(Dispersal~magnitude~(m^2~day^-1)),expand=c(0,0),
                     labels=labs_kh,breaks=paste("E",1:length(k_h),sep=""))+
    theme_bw()+ 
    lat_div+
    coord_cartesian(clip="off",ylim = c(1,159))+
    geom_text(data=data_text, aes(x=0.5, y=164, label=label), fontface="bold", size=4)+
    facet_grid(~factor.m,labeller=label_parsed)

  ggsave(paste(figure_home,"div_disp-",subtitle,".png",sep=""),dpi=600,plot.gdiv,width=6.5,height=4.5, units = "in")
  
  
  
  ### Niche parameters ----
  vals.all<-data.frame(do.call("rbind",lapply(spatial_analysis,function(y) {y$vals})))
  for (i in 1:n) {
    vals.all$t.rat[vals.all$name==i]<-(vals.all$tobs[vals.all$name==i]-fund_topt[i])/fund_topt[i]
    vals.all$w.rat[vals.all$name==i]<-(vals.all$wobs[vals.all$name==i]-fund_widths[i])/fund_widths[i]
  }
  vals.all$ogmag<-as.numeric(str_extract(vals.all$mag,"\\d"))
  vals.all$mag<- as.numeric(str_extract(vals.all$mag,"\\d"))
  for(f in 1:length(mort)) {
    jitter<-seq(-0.2,0.2,length_out = length(mort))[f]
    vals.all$mag[vals.all$m==paste0("exp",f)]<-vals.all$mag[vals.all$m==paste0("exp",f)]+jitter
  }
  vals.all$factor.m<-factor(vals.all$m,levels= levels(vals.all$m)[order(mort)])
  vals.all<-vals.all%>%filter(mag<4.5)
  
  n.par<-theme(panel.grid = element_blank(),
               text = element_text(color="black",size=10),
               panel.border = element_rect(color="black"),
               axis.title = element_text(size=10),
               axis.text = element_text(size=9,color="black"),
               axis.text.x = element_text(angle = 25,vjust=1,hjust=1),
               plot.margin = unit(c(1,0.5,0.5,0.5),"lines"))
  
  spatial.glm.t<-glm(data=vals.all,formula=t.rat~bs(ogmag)*factor.m)
  spatial.glm.w<-glm(data=vals.all,formula=w.rat~bs(ogmag)*factor.m)
   summary(spatial.glm.w)

  plot.a<-ggplot(data=vals.all,aes(x=jitter(mag),y=t.rat, color=factor.m))+
    geom_point(alpha=0.3)+
    geom_line(aes(x=mag,y=predict(spatial.glm.t)))+
    geom_boxplot(aes(group=mag, fill=factor.m),color="black", size=0.2, outlier.color = NA,show.legend = F)+
    
    scale_y_continuous(name=bquote(delta[T]))+
    geom_hline(yintercept=0,linetype="dashed")+
    scale_x_continuous(name=bquote(Dispersal~magnitude~(m^2~day^-1)),expand=c(0.005,0.005),
                       labels=labs_kh,breaks=unique(round(vals.all$mag)))+
    scale_color_manual(name=bquote(italic(m)~(day^-1)),
                       values=niche_colors,
                       labels=as.numeric(mort)[order(mort)])+
    scale_fill_manual(name=bquote(Mortality~(day^-1)),
                       values=niche_colors,
                       labels=as.numeric(mort)[order(mort)])+
    theme_bw()+n.par+theme(legend.position = "none",
                           axis.title.y = element_text(size=13)) #+ntpar
  
  plot.b<-ggplot(data=vals.all,aes(x=jitter(mag),y=w.rat, color=factor.m))+
    geom_point(alpha=0.3)+
    geom_line(aes(x=mag,y=predict(spatial.glm.w)))+
    geom_boxplot(aes(group=mag, fill=factor.m),size=0.2,color="black", outlier.color = NA,show.legend = F)+
    
  #geom_smooth(method="glm",aes(x=mag,y=w.rat,col=factor.m,group=factor.m),se=F)+
  geom_hline(yintercept=0,linetype="dashed")+
    scale_y_continuous(name=bquote(delta[W]))+
    scale_x_continuous(name=bquote(Dispersal~magnitude~(m^2~day^-1)),expand=c(0.005,0.005),
                       labels=labs_kh,breaks=unique(round(vals.all$mag)))+
    scale_color_manual(name=bquote(Mortality~(day^-1)),
                       values=niche_colors,
                       labels=as.numeric(mort)[order(mort)])+
    scale_fill_manual(name=bquote(Mortality~(day^-1)),
                       values=niche_colors,
                       labels=as.numeric(mort)[order(mort)])+
    theme_bw()+n.par+theme(legend.position = c(0.2,0.88), 
                           axis.title.y = element_text(size=13),
                           legend.background = element_rect("transparent"),
                           legend.spacing.y = unit(0,"lines"),
                           legend.title = element_text(size=9),
                           legend.key.height = unit(0.1, "in"))+
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
  niche.disp<-ggarrange(plot.b,plot.a, labels = c("a.","b.")) #,label.y=1,label.x = 0.05,common.legend = T,legend = "right")
  ggsave(paste(figure_home,"niche_disp-boxplot",subtitle,".png",sep=""),dpi=600,niche.disp,width=6.5,height=3.5)
  


  
  ### Summary Values ----
  alpha.all%>%
    group_by(m,mag)%>%
    summarise(mean=mean(value), min=min(value),max=max(value))
  summary(aov(value~m*factor(mag), data=alpha.all))

  vals.all%>%
    group_by(m,mag)%>%
    summarise(topt=mean(tobs),w=mean(wobs),tdiff=mean(t.rat),wdiff=mean(w.rat))
  summary(aov(w.rat~factor.m*factor(ogmag), data=vals.all))
  summary(aov(t.rat~factor.m*factor(ogmag), data=vals.all))
  
#### E3: Temporal Mass Effects Figures -----
  ### Proof of Concept Plot ----
  only_t<-str_detect(file_list, "T1_D0_E0_exp1")
  load(paste(files_home,file_list[only_t],sep=""))

  div.list.temp<-temporal_analysis[[1]][[2]]$species

  
  l<-boxes_to_plot
  p1<-do.call("rbind",lapply(1:length(l),function(c) {
    if(length(unique(unlist(div.list.temp[[l[c]]])))==1) {
      df<-data.frame(box=l[c],time_vec=as.vector(time_vec[lastyear]/365),temp=as.vector(TEMP[lastyear,l[c]]),P[lastyear,unlist(div.list.temp[[l[c]]]),l[c]])
      colnames(df)[4]<-paste0("Species.",unlist(div.list.temp[[l[c]]])) 
    }else{
      df<-data.frame(box=l[c],time_vec=as.vector(time_vec[lastyear]/365),temp=as.vector(TEMP[lastyear,l[c]]),P[lastyear,unlist(div.list.temp[[l[c]]]),l[c]])
    }
    df%>%pivot_longer(cols = !c("box","time_vec","temp"))
  }))
  plot.df<-p1
  plot.df$name<-factor(plot.df$name,levels=paste("Species",1:n,sep="."))
  plot.df$lt<-"np"
  plot.df$lt[plot.df$name %in% sp.names.none & plot.df$value>1E-3]="p"
  
  labels.temp<-paste(letters[seq(1,length(l)*2,by=2)],".",sep="")
  names(labels.temp)<-l
  labels.ab<-paste(letters[seq(2,(length(l)*2),by=2)],".",sep="")
  names(labels.ab)<-l
  plot.a<-ggplot(plot.df)+
    geom_line(aes(x=time_vec,y=temp))+
    scale_y_continuous(name=bquote(Temperature~(degree*C)))+
    scale_x_continuous(name=bquote(time_vec~(year)),expand=c(0,0),breaks=quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),labels=signif(quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),3))+
    facet_wrap(~box,ncol = 1,labeller = labeller(box=labels.temp))+
    theme_bw()+poc.par
  sp.names<-unique(plot.df$name)[order(unique(plot.df$name))]
  sp.use.plot<-gsub("Species.","",sp.names)
  
  sp.names.temporal<-sp.names
  
  sp.use.plot.2<-z_vals[as.numeric(sp.use.plot)][order(z_vals[as.numeric(sp.use.plot)])]
  names(sp.use.plot.2)<-sp.names[order(as.numeric(sp.use.plot))]
  plot.b<-ggplot(plot.df)+
    geom_line(aes(x=time_vec,y=log10(value),col=name,group=name,linetype=lt))+
    scale_color_manual(name=bquote(T[opt]~(degree*C)),
                       labels=sp.use.plot.2,values=colors_P(length(sp.use.plot.2)))+
    scale_y_continuous(name=bquote(log[10]~Biomass~(mmol~P~m^-3)),limits=c(-4,1.5),expand=c(0,0))+
    scale_x_continuous(name=bquote(time_vec~(year)),expand=c(0,0),breaks=quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),labels=signif(quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),3))+
    facet_wrap(~box,ncol = 1,labeller = labeller(box=labels.ab))+
    scale_linetype_manual(values=c("np"="solid","p"="dotdash"),guide=F)+
    theme_bw()+poc.par
  
  
  exp_temporal=list(a=plot.a,b=plot.b)
  save(exp_temporal,file=paste0(proof_files,"exp-temp_plot-",subtitle,".Rdata"))
  
  
  load(paste0(files_home,"exp-temp_plot-",subtitle,".Rdata"))
  proof.temp<-ggarrange(exp_temporal$a,exp_temporal$b,ncol=2,common.legend = T,legend="right")
  ggsave(paste(figure_home,"proof_temp.svg",sep=""),dpi=600,proof.temp,width=4,height=6)
  
  
  
  
  ### Diversity Plot -----
  ### 
  #TEMP<-TEMP.in
  gam.div<-as.data.frame(do.call("rbind",lapply(temporal_analysis,function(y) {y$div$gamma})))
  colnames(gam.div)<-apply(TEMP.in,2,function(x) {diff(range(x))})
  rownames(gam.div)<-paste("m",1:length(mort),sep="_")
  gam.div<-gam.div%>%
    mutate(m=rownames(.))%>%
    pivot_longer(colnames(gam.div)[!colnames(gam.div)=="m"])
  
  gam.div$factor.m<-factor(gam.div$m,levels= unique(gam.div$m)[order(mort)])
  
  alpha.div<-as.data.frame(do.call("rbind",lapply(temporal_analysis,function(y) {y$div$alpha})))
  colnames(alpha.div)<-apply(TEMP.in,2,function(x) {diff(range(x))})
  rownames(alpha.div)<-paste("m",1:length(mort),sep="_")
  alpha.div<-alpha.div%>%
    mutate(m=rownames(.))%>%
    pivot_longer(colnames(alpha.div)[!colnames(alpha.div)=="m"])
  
  alpha.div$factor.m<-factor(alpha.div$m,levels= unique(alpha.div$m)[order(mort)])
  
  tamp.par<-theme(panel.grid = element_blank(),
               text = element_text(color="black",size=10),
               panel.border = element_rect(color="black"),
               axis.title = element_text(size=10),
               axis.text = element_text(size=9,color="black"),
               plot.margin = unit(c(1,0.5,0.5,0.5),"lines"))
  
  temp.glm.a<-gam(data=alpha.div,formula=value~bs(as.numeric(name))*m)
  temp.glm.g<-gam(data=gam.div,formula=value~bs(as.numeric(name))*m)
  
  plot.a<-ggplot(gam.div, aes(x=as.numeric(name),y=value,col=factor.m, group=factor.m))+
    geom_point()+
    geom_line(aes(y=predict.Gam(temp.glm.g),color=m))+
    scale_y_continuous(name=bquote(S[T]))+
    scale_x_continuous(name=bquote(Temperature~amplitude~(degree*C)))+
    scale_color_manual(name=bquote(Mortality~(day^-1)),
                       labels=mort[order(mort)],
                       values=niche_colors)+
    coord_cartesian(expand=0,xlim=c(0,10.1), ylim=c(0.8,5.2))+
    theme_bw()+tamp.par+theme(legend.position = "none")
  
  plot.b<-ggplot(alpha.div, aes(x=as.numeric(name),y=value,col=factor.m, group=factor.m))+
    geom_point()+
    geom_line(aes(y=predict.Gam(temp.glm.a),color=m))+
    scale_y_continuous(name=bquote(bar(S)))+
    scale_x_continuous(name=bquote(Temperature~amplitude~(degree*C)))+
    scale_color_manual(name=bquote(Mortality~(day^-1)),
                       labels=mort[order(mort)],
                       values=niche_colors)+
    coord_cartesian(expand=0,xlim=c(0,10.1), ylim=c(0.8,5.2))+
    theme_bw()+tamp.par+theme(legend.position = c(0.2,0.88), 
                           legend.background = element_rect("transparent"),
                           legend.spacing.y = unit(0,"lines"),
                           legend.title = element_text(size=9),
                           legend.key.height = unit(0.1, "in"))
  
  div.temp<-ggarrange(plot.b, plot.a,
                      labels = c("a.","b.")) #,label.x = 0.1,label.y = 1.01,common.legend = T,legend="right")
  ggsave(paste(figure_home,"div_temp-",subtitle,".png",sep=""),dpi=600,div.temp,width=6.5,height=3.5)
  
  #### Try something new
  # gam.div$div="gamma"
  # alpha.div$div="alpha"
  # div.disp<-rbind(gam.div,alpha.div)
  # div.combined.temp<-ggplot(subset(div.disp,m=="m_1"))+
  #   geom_point(aes(x=as.numeric(name),y=value,col=div))+
  #   stat_smooth(aes(x=as.numeric(name),y=value,col=div,group=div),method = "loess",se=F)+
  #   scale_y_continuous(name=bquote(Number~of~species))+
  #   scale_x_continuous(name=bquote(Temperature~amplitude~(degree*C)))+
  #   scale_color_manual(name="",
  #                      labels=c(bquote(bar(S)),bquote(S[T])),
  #                      values=amp_div_colors)+
  #   coord_cartesian(expand=0,xlim=c(0,10.1), ylim=c(0.8,5.2))+
  #   theme_bw()+
  #   tamp.par
  # 
  # ggsave(paste(figure_home,"div_temp_combo-",subtitle,".png",sep=""),dpi=600,div.combined.temp,width=5,height=5)
  
  
  
  ### Niche parameters ----
  vals.all.temp<-data.frame(do.call("rbind",lapply(temporal_analysis,function(y) {y$vals})))
  for (i in 1:n) {
    vals.all.temp$t.rat[vals.all.temp$name==i]<-(vals.all.temp$tobs[vals.all.temp$name==i]-fund_topt[i])/fund_topt[i]
    vals.all.temp$w.rat[vals.all.temp$name==i]<-(vals.all.temp$wobs[vals.all.temp$name==i]-fund_widths[i])/fund_widths[i]
    vals.all.temp$topt[vals.all.temp$name==i]<-fund_topt[i]
  }
  vals.all.temp$factor.m<-factor(vals.all.temp$m,levels= levels(vals.all.temp$m)[order(mort)])
  
  glmt<-glm(formula = t.rat~m+tdiff,data=vals.all.temp)
  glmw<-glm(formula = w.rat~m+tdiff,data=vals.all.temp)
  # summary(glht(glmw, mcp(m="Tukey")))
  # summary(glht(glmt, mcp(m="Tukey")))
  
  
  plot.a<-ggplot(data=vals.all.temp, aes(x=tdiff,y=t.rat,col=factor.m))+
    geom_point()+
    geom_line(aes(y=predict(glmt)))+
    #geom_smooth(method="loess", formula= y~x,aes(x=tdiff,y=t.rat,col=factor.m,group=factor.m),se=F)+
    scale_y_continuous(name=bquote(delta[T]))+
    geom_hline(yintercept=0,linetype="dashed")+
    scale_x_continuous(name=bquote(Temperature~amplitude~(degree*C)))+
    scale_color_manual(name=bquote(italic(Mortality)~(day^-1)),
                       values=niche_colors,
                       labels=as.numeric(mort)[order(mort)])+
    theme_bw()+
    tamp.par+theme(legend.position = "none",axis.title.y = element_text(size=13))
  
  
  plot.b<-ggplot(data=vals.all.temp,aes(x=tdiff,y=w.rat,col=factor.m))+
    geom_point()+
    geom_line(aes(y=predict(glmw)))+
    #geom_smooth(method="loess", formula= y~x,aes(x=tdiff,y=w.rat,col=factor.m,group=factor.m),se=F)+
    geom_hline(yintercept=0,linetype="dashed")+
    scale_y_continuous(name=bquote(delta[W]))+
    scale_x_continuous(name=bquote(Temperature~amplitude~(degree*C)))+
    scale_color_manual(name=bquote(Mortality~(day^-1)),
                       values=niche_colors,
                       labels=as.numeric(mort)[order(mort)])+
    theme_bw()+
    tamp.par+theme(legend.position = c(0.2,0.88), 
           legend.background = element_rect("transparent"),
           legend.spacing.y = unit(0,"lines"),
           axis.title.y = element_text(size=13),
           legend.title = element_text(size=9),
           legend.key.height = unit(0.1, "in"))
  
  niche.temp<-ggarrange(plot.b,plot.a,labels = paste0(letters[1:2],".")) #,label.y=1,label.x = 0.05,common.legend = T,legend = "right")
  
  ggsave(paste(figure_home,"niche_temp-",subtitle,".png",sep=""),dpi=600,niche.temp,width=6.5,height=3.5)
  
  
  
  ### Summary Values ----
  gam.div%>%
    group_by(m)%>%
    summarise(mean=mean(value),range(value))
  alpha.div%>%
    group_by(m)%>%
    summarise(mean=mean(value),range(value))
  vals.all.temp%>%
    group_by(m)%>%
    summarise(topt=mean(tobs),w=mean(wobs),tdiff=mean(t.rat),wdiff=max(w.rat), min=max(wobs))
  
#### E4: Combo Effects Figures -----
  ### Proof of Concept Plot ----
  only_t<-str_detect(file_list, "T1_D1_E1_exp1")
  load(paste(files_home,file_list[only_t],sep=""))

  all_output<-combined_analysis[[1]][[2]]$species

  
  l<-boxes_to_plot
  p1<-do.call("rbind",lapply(1:length(l),function(c) {
    if(length(unique(unlist(all_output[[l[c]]])))==1) {
      df<-data.frame(box=l[c],time_vec=as.vector(time_vec[lastyear]/365),temp=as.vector(TEMP[lastyear,l[c]]),P[lastyear,unlist(all_output[[l[c]]]),l[c]])
      colnames(df)[4]<-paste0("Species.",unlist(all_output[[l[c]]]))
    }else{
      df<-data.frame(box=l[c],time_vec=as.vector(time_vec[lastyear]/365),temp=as.vector(TEMP[lastyear,l[c]]),P[lastyear,unlist(all_output[[l[c]]]),l[c]])
    }
    df%>%pivot_longer(cols = !c("box","time_vec","temp"))
  }))
  plot.df<-p1
  plot.df$name<-factor(plot.df$name,levels=paste("Species",1:n,sep="."))
  plot.df$lt<-"np"
  plot.df$lt[plot.df$name %in% sp.names.none]="p"

  labels.temp<-paste(letters[seq(1,length(l)*2,by=2)],".",sep="")
  names(labels.temp)<-l
  labels.ab<-paste(letters[seq(2,(length(l)*2),by=2)],".",sep="")
  names(labels.ab)<-l
  plot.a<-ggplot(plot.df)+
    geom_line(aes(x=time_vec,y=temp))+
    scale_y_continuous(name=bquote(Temperature~(degree*C)))+
    scale_x_continuous(name=bquote(time_vec~(year)),expand=c(0,0),breaks=quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),labels=signif(quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),3))+
    facet_wrap(~box,ncol = 1,labeller = labeller(box=labels.temp))+
    theme_bw()+poc.par
  sp.names<-unique(plot.df$name)[order(unique(plot.df$name))]
  sp.use.plot<-gsub("Species.","",sp.names)
  sp.use.plot.2<-z_vals[as.numeric(sp.use.plot)][order(z_vals[as.numeric(sp.use.plot)])]
  names(sp.use.plot.2)<-sp.names[order(as.numeric(sp.use.plot))]
  plot.b<-ggplot(plot.df)+
    geom_line(aes(x=time_vec,y=log10(value),col=name,group=name,linetype=lt))+
    scale_color_manual(name=bquote(T[opt]~(degree*C)),
                       labels=sp.use.plot.2,values=colors_P(length(sp.use.plot.2)))+
    scale_linetype_manual(values=c("np"="solid","p"="dotdash"),guide=F)+
    scale_y_continuous(name=bquote(log[10]~Biomass~(mmol~P~m^-3)),limits=c(-4,1.5),expand=c(0,0))+
    scale_x_continuous(name=bquote(time_vec~(year)),expand=c(0,0),breaks=quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),labels=signif(quantile(as.vector(time_vec[lastyear]/365),c(0,0.5,1)),3))+
    facet_wrap(~box,ncol = 1,labeller = labeller(box=labels.ab))+
    theme_bw()+poc.par
  
  exp_combo=list(a=plot.a,b=plot.b)
  save(exp_combo,file=paste0(proof_files,"exp-combo_plot-",subtitle,".Rdata"))
  
  
  
  load(paste0(files_home,"exp-combo_plot-",subtitle,".Rdata"))
  proof.combo<-ggarrange(exp_combo$a,exp_combo$b,ncol=2,common.legend = T,legend="right")
  ggsave(paste(figure_home,"proof_combo.svg",sep=""),dpi=600,proof.combo,width=4.5,height=6)
  
  ###Just testing something
  # testing.subdf<-subset(plot.df,name=="Species.24" & box==l[2])
  # cov(testing.subdf$value,mu.t[lastyear,27,l[2]])
  
  
  ### Diversity plot ----
  ### 
    load("~/Desktop/combined_analysis-testing_30w.Rdata")
  gam.div.combo<-as.data.frame(do.call("rbind",lapply(combined_analysis,function(x) {c(x$div$gamma,x$div$ids)})))
  colnames(gam.div.combo)<-c(paste("Lat.",1:bx,sep=""),"T","disp","mag","m") #c(apply(TEMP,2,function(x) {diff(range(x))}),"T","disp","mag","m")
  gam.div.combo<-gam.div.combo%>%
    pivot_longer(colnames(gam.div.combo)[!colnames(gam.div.combo)%in%c("m","T","disp","mag")])
  gam.div.combo$factor.m<-factor(gam.div.combo$m,levels=unique(gam.div.combo$m),labels=mort_labs) #paste(c("bold(a.)~","bold(c.)~"),mort_labs))
  gam.div.combo$div<-"S[T]"
  alpha.div.combo<-as.data.frame(do.call("rbind",lapply(combined_analysis,function(x) {c(x$div$alpha,x$div$ids)})))
  colnames(alpha.div.combo)<-c(paste("Lat.",1:bx,sep=""),"T","disp","mag","m") #c(apply(TEMP,2,function(x) {diff(range(x))}),"T","disp","mag","m")
  alpha.div.combo<-alpha.div.combo%>%
    pivot_longer(colnames(alpha.div.combo)[!colnames(alpha.div.combo)%in%c("m","T","disp","mag")])
  alpha.div.combo$factor.m<-factor(alpha.div.combo$m,levels=unique(alpha.div.combo$m),labels=mort_labs)
  alpha.div.combo$div="bar(S)"
  
  total.div<-rbind(alpha.div.combo,gam.div.combo)
  total.div$name<-factor(total.div$name,levels=unique(total.div$name))
  total.div$factor.m<-factor(total.div$factor.m,levels=levels(total.div$factor.m)[order(levels(total.div$factor.m))])
  total.div<-total.div%>%filter(!mag=="E5")
  
  data_text <- data.frame(
    factor.m = rep(mort, 2),
    div = rep(unique(total.div$div), each=3),
    label = paste0(letters[1:6],".")
  )
  data_text$factor.m<-factor(data_text$factor.m,levels=unique(data_text$factor.m),labels=mort_labs)
  #data_text$factor.m<-factor(data_text$factor.m,levels=levels(data_text$factor.m)[order(mort)])
  
  total.div$IAV<-rep(apply(TEMP.in,2,function(x) {diff(range(x))}),2*length(k_h)*length(mort))
  
  
  combo.glms.alpha<-gam(formula=as.numeric(as.vector(value))~mag*IAV+mag*m,data=subset(total.div,div=="bar(S)"))
  combo.glms.gamma<-gam(formula=as.numeric(as.vector(value))~mag*IAV+mag*m,data=subset(total.div,div=="S[T]"))
  
  summary(combo.glms.alpha)
  summary(combo.glms.gamma)
  
  summary(glht(combo.glms.alpha, mcp(m="Tukey")))
  summary(glht(combo.glms.gamma, mcp(m="Tukey")))
  
  
  
  combo.heat<-ggplot(total.div)+
    geom_tile(aes(x=mag,y=name,fill=as.numeric(as.character(value))))+
    scale_fill_gradientn(name="# Species",colours=brewer.pal(9,"Greens")[2:9],na.value = "grey",#limits=c(0.01,NA),
                         breaks=pretty(seq(1,max(as.numeric(as.character(alpha.all$value))),length_out=5)))+
    scale_y_discrete(name=bquote(Latitude),expand=c(0,0),
                     labels=lats_plot,
                     breaks=paste("Lat",which(full_latitude%in%lats_plot),sep="."))+
    scale_x_discrete(name=bquote(Dispersal~magnitude~(m^2~day^-1)),expand=c(0,0),
                     labels=labs_kh,breaks=paste("E",1:length(k_h),sep=""))+
    theme_bw()+ 
    lat_div+
    theme(panel.spacing.y = unit(1,"lines"))+
    coord_cartesian(clip="off",ylim = c(1,159))+
    geom_text(data=data_text, aes(x=0.5, y=170, label=label), fontface="bold", size=4)+
    facet_grid(div~factor.m,labeller=label_parsed)
    
  
  ggsave(paste(figure_home,"combo_div_heat-",subtitle,".png",sep=""),dpi=600,combo.heat,width=6.5,height=4.5)
  
  
  ### Niche parameters ----
  all.values<-data.frame(do.call("rbind",lapply(combined_analysis,function(y) {y$vals})))
  local.values<-data.frame(do.call("rbind",lapply(combined_analysis,function(y) {y$vals.box})))
  for (i in 1:n) {
    all.values$t.rat[all.values$name==i]<-(all.values$tobs[all.values$name==i]-fund_topt[i])/fund_topt[i]
    all.values$w.rat[all.values$name==i]<-(all.values$wobs[all.values$name==i]-fund_widths[i])/fund_widths[i]
    local.values$t.rat[local.values$name==i]<-(local.values$tobs[local.values$name==i]-fund_topt[i])/fund_topt[i]
    local.values$w.rat[local.values$name==i]<-(local.values$wobs[local.values$name==i]-fund_widths[i])/fund_widths[i]
  }
  local.values$m<-factor(local.values$m,levels=pattern_search,labels=mort)
  all.values$m<-factor(all.values$m,levels=pattern_search,labels=mort)
  local.values$factor.m<-factor(local.values$m,levels=mort,labels=mort_labs)
  all.values$factor.m<-factor(all.values$m,levels=mort,labels=mort_labs)
  #### Jitterin g
  #### 
  all.values$mag<- as.numeric(str_extract(all.values$mag,"\\d"))
  for(f in 1:length(mort)) {
    jitter<-seq(-0.2,0.2,length_out = length(mort))[f]
    all.values$mag[all.values$m==mort[order(mort)][f]]<-all.values$mag[all.values$m==mort[order(mort)][f]]+jitter
  }
  
  all.values$m<-factor(all.values$m,levels=levels(all.values$m)[order(mort)])
  all.values<-all.values%>%filter(mag<4.5)
  
  combo.glms.w<-gam(formula=w.rat~mag*tdiff+m,data=all.values)
  combo.glms.t<-gam(formula=t.rat~mag*tdiff+m,data=all.values)
  combo.glms.wonly<-glm(formula=w.rat~mag*m,data=all.values)
  combo.glms.tonly<-glm(formula=t.rat~mag*m,data=all.values)
  
  summary(combo.glms.wonly)
  summary(combo.glms.t)
  
  summary(glht(combo.glms.w, mcp(m="Tukey")))
  summary(glht(combo.glms.t, mcp(m="Tukey")))
  
  
  t.disp<-ggplot(data=all.values,aes(x=mag,y=w.rat,fill=tdiff,shape=m))+
    geom_point(color="white",size=4)+
     geom_line(aes(y=predict(combo.glms.wonly),linetype=m))+
    #geom_smooth(method="lm", formula= y~x,aes(x=round(mag),y=w.rat,linetype=m),se=F,col="black",size=0.8)+
    geom_hline(yintercept=0,linetype="dashed")+
    scale_fill_gradientn(name=bquote(Delta*T),colors=rev(brewer.pal(9,"Spectral")))+
    scale_y_continuous(name=bquote(delta[W]))+
    scale_x_continuous(name=bquote(Dispersal~magnitude~(m^2~day^-1)),expand=c(0.1,0.1),
                       labels=labs_kh,breaks=unique(round(all.values$mag)))+
    scale_shape_manual(name=bquote(Mortality~(day^-1)),
                       labels=mort[order(mort)],values=c(21,22,24),
                       guide=guide_legend(override.aes = list(color="black")))+
    scale_linetype_manual(name=bquote(Mortality~(day^-1)),
                            labels=mort[order(mort)],values=c(1,3,6))+
    # scale_color_manual(name=bquote(Mortality~(day^-1)),
    #                    values=niche_colors,
    #                    labels=mort[order(mort)])+
    theme_bw()+n.par+theme(axis.title.y = element_text(size=13))
  
  
  #ggsave("~/Documents/combined_width.png",dpi=600,t.disp,width=5,height=4.5)
  t.disp.t<-ggplot(data=all.values,aes(x=mag,y=t.rat,fill=tdiff,shape=m))+
    geom_point(color="white",size=4)+
    geom_line(aes(y=predict(combo.glms.tonly),linetype=m))+
    #geom_smooth(method="lm", formula= y~x,aes(x=round(mag),y=w.rat,linetype=m),se=F,col="black",size=0.8)+
    geom_hline(yintercept=0,linetype="dashed")+
    scale_fill_gradientn(name=bquote(Delta*T),colors=rev(brewer.pal(9,"Spectral")))+
    scale_y_continuous(name=bquote(delta[T]))+
    scale_x_continuous(name=bquote(Dispersal~magnitude~(m^2~day^-1)),expand=c(0.1,0.1),
                       labels=labs_kh,breaks=unique(round(all.values$mag)))+
    scale_shape_manual(name=bquote(Mortality~(day^-1)),
                       labels=mort[order(mort)],values=c(21,22,24),
                       guide=guide_legend(override.aes = list(color="black")))+
    scale_linetype_manual(name=bquote(Mortality~(day^-1)),
                            labels=mort[order(mort)],values = c(1,3,6))+
    # scale_color_manual(name=bquote(Mortality~(day^-1)),
    #                    values=niche_colors,
    #                    labels=mort[order(mort)])+
    theme_bw()+n.par+theme(axis.title.y = element_text(size=13))
  
  t.disp.combo<-ggarrange(t.disp,t.disp.t,common.legend = T,legend = "right",labels = c("a.","b."),label.y=1,label.x = 0.01)+
    theme(plot.background = element_rect(fill="white", color = NA))

  ggsave(paste(figure_home,"combo_niche-",subtitle,"lines.png",sep=""),dpi=600,t.disp.combo,width=7.5,height=3.5)
  
   
  ### Summary Values ----
  total.div%>%
    group_by(factor.m,div,mag)%>%
    summarise(mean=mean(as.numeric(as.vector((value)))),min=max(as.numeric(as.vector((value)))))%>%
    as.data.frame()
  all.values%>%
    group_by(m,mag)%>%
    summarise(topt=mean(tobs),w=mean(wobs),tdiff=mean(t.rat),wdiff=mean(w.rat))
  
#### Comparing net rates ------
  
  ### Growth
  
  only_t<-str_detect(file_list, "T1_D1_E3_exp1")
  load(paste(files_home,file_list[only_t],sep=""))
  #load(paste(files_home,"combined_analysis-",subtitle,".Rdata",sep=""))
  
  div.spa<-combined_analysis[[3]]$div$species
  
  taumatrix0=matrix(0,nrow=bx,ncol=bx)
  ######Name the rows and columns amu.ter the box numbers
  colnames(taumatrix0)=full_latitude
  rownames(taumatrix0)=full_latitude
  ######Take the difference between the rows and columns for the difference factor
  diffmatrix=abs(outer(as.numeric(colnames(taumatrix0)),as.numeric(rownames(taumatrix0)),"-"))
  colnames(diffmatrix)=full_latitude
  rownames(diffmatrix)=full_latitude
  ####Give realistic differences (1 degree full_latitude= 111km convert to m)
  diffmatrix<-diffmatrix*(111*1000)
  ########Create the tau matrix 
  
  taumatrix2=k_h[3]/(diffmatrix^2)
  taumatrix2[which(!is.finite(taumatrix2))]=0
  
t.mortality<-(0.81*exp(0.0631*TEMP))

png(paste0(figure_home,"testing_nets_3-tvar.png"),width = 12, height=10, units="in", res = 300)
species=10 #20 #33
boxes<-1:bx #which(unlist(lapply(div.spa,function(u) {species%in%u})))
growth<-(mu.t[lastyear,species,boxes]*P[lastyear,species,boxes]*(R[lastyear,species,boxes]/(R[lastyear,species,boxes]+k)))
const.mortality<-0.05*P[lastyear,species,boxes]
decay<-const.mortality*t.mortality[lastyear,boxes]

gmdiff<-apply(growth-decay,2,function(y) {as.vector(by(y, ceiling(lastyear / oneyr_iter), function(z) {mean(z,na.rm=T)}))})[2,]
Pdiff<-apply(P[lastyear,species,boxes],2,function(y) {as.vector(by(y, ceiling(lastyear / oneyr_iter), function(z) {mean(z,na.rm=T)}))})[2,]

P2=array(apply(P[lastyear,species,boxes],1,function(x) {rep(x,each=length(boxes))}),dim=c(length(boxes),length(boxes),length(lastyear))) 
tausum<-t(apply(P2,3,function(x) {apply(x*taumatrix2[boxes,boxes],1,function(y) {sum(y,na.rm=T)})}))
tausub<-t(apply(P2,3,function(x) {apply(x*taumatrix2[boxes,boxes],2,function(y) {sum(y,na.rm=T)})}))

taudiff<-apply(tausum-tausub,2,function(y) {as.vector(by(y, ceiling(lastyear / oneyr_iter), function(z) {mean(z,na.rm=T)}))})[2,]
#taudiff<-rep(0,bx)

par(mar=c(4, 8, 2, 6) + 0.1,mfrow=c(3,1))
plot(log10(Pdiff), axes=F, xlab="",ylab="",type="l",col="black", main=bquote(T[opt]==.(z_vals[species])*degree*C), ylim=c(-6,1))
lines(rep(-4,length(Pdiff)),lty=2, col="black")
axis(4,col="black",lwd=2)
mtext(4,text="Biomass", line=2)
par(new=T)
plot(taudiff, axes=F, xlab="",ylab="",type="l",col="blue", main="")
lines(rep(0,length(taudiff)),lty=2, col="blue")
axis(2,col="blue",lwd=2)
mtext(2,text="Net Dispersal",line=2)
par(new=T)
plot(gmdiff, axes=F, xlab="",ylab="",type="l",col="red", main="")
lines(rep(0,length(gmdiff)),lty=2, col="red")
axis(2,col="red",lwd=2, line=3.5)
mtext(2,text="mu-m",line=5.5)
axis(1,at=1:bx,labels = full_latitude)
mtext("Latitude",side=1,col="black",line=2)


species=20 #20 #33
boxes<-1:bx #which(unlist(lapply(div.spa,function(u) {species%in%u})))
growth<-(mu.t[lastyear,species,boxes]*P[lastyear,species,boxes]*(R[lastyear,species,boxes]/(R[lastyear,species,boxes]+k)))
const.mortality<-0.05*P[lastyear,species,boxes]
decay<-const.mortality*t.mortality[lastyear,boxes]

gmdiff<-apply(growth-decay,2,function(y) {as.vector(by(y, ceiling(lastyear / oneyr_iter), function(z) {mean(z,na.rm=T)}))})[2,]
Pdiff<-apply(P[lastyear,species,boxes],2,function(y) {as.vector(by(y, ceiling(lastyear / oneyr_iter), function(z) {mean(z,na.rm=T)}))})[2,]

P2=array(apply(P[lastyear,species,boxes],1,function(x) {rep(x,each=length(boxes))}),dim=c(length(boxes),length(boxes),length(lastyear))) 
tausum<-t(apply(P2,3,function(x) {apply(x*taumatrix2[boxes,boxes],1,function(y) {sum(y,na.rm=T)})}))
tausub<-t(apply(P2,3,function(x) {apply(x*taumatrix2[boxes,boxes],2,function(y) {sum(y,na.rm=T)})}))

taudiff<-apply(tausum-tausub,2,function(y) {as.vector(by(y, ceiling(lastyear / oneyr_iter), function(z) {mean(z,na.rm=T)}))})[2,]
#taudiff<-rep(0,bx)

#par(mar=c(5, 8, 4, 6) + 0.1)
plot(log10(Pdiff), axes=F, xlab="",ylab="",type="l",col="black", main=bquote(T[opt]==.(z_vals[species])*degree*C), ylim=c(-6,1))
lines(rep(-4,length(Pdiff)),lty=2, col="black")
axis(4,col="black",lwd=2)
mtext(4,text="Biomass", line=2)
par(new=T)
plot(taudiff, axes=F, xlab="",ylab="",type="l",col="blue", main="")
lines(rep(0,length(taudiff)),lty=2, col="blue")
axis(2,col="blue",lwd=2)
mtext(2,text="Net Dispersal",line=2)
par(new=T)
plot(gmdiff, axes=F, xlab="",ylab="",type="l",col="red", main="")
lines(rep(0,length(gmdiff)),lty=2, col="red")
axis(2,col="red",lwd=2, line=3.5)
mtext(2,text="mu-m",line=5.5)
axis(1,at=1:bx,labels = full_latitude)
mtext("Latitude",side=1,col="black",line=2)

species=33 #20 #33
boxes<-1:bx #which(unlist(lapply(div.spa,function(u) {species%in%u})))
growth<-(mu.t[lastyear,species,boxes]*P[lastyear,species,boxes]*(R[lastyear,species,boxes]/(R[lastyear,species,boxes]+k)))
const.mortality<-0.05*P[lastyear,species,boxes]
decay<-const.mortality*t.mortality[lastyear,boxes]

gmdiff<-apply(growth-decay,2,function(y) {as.vector(by(y, ceiling(lastyear / oneyr_iter), function(z) {mean(z,na.rm=T)}))})[2,]
Pdiff<-apply(P[lastyear,species,boxes],2,function(y) {as.vector(by(y, ceiling(lastyear / oneyr_iter), function(z) {mean(z,na.rm=T)}))})[2,]

P2=array(apply(P[lastyear,species,boxes],1,function(x) {rep(x,each=length(boxes))}),dim=c(length(boxes),length(boxes),length(lastyear))) 
tausum<-t(apply(P2,3,function(x) {apply(x*taumatrix2[boxes,boxes],1,function(y) {sum(y,na.rm=T)})}))
tausub<-t(apply(P2,3,function(x) {apply(x*taumatrix2[boxes,boxes],2,function(y) {sum(y,na.rm=T)})}))

taudiff<-apply(tausum-tausub,2,function(y) {as.vector(by(y, ceiling(lastyear / oneyr_iter), function(z) {mean(z,na.rm=T)}))})[2,]
#taudiff<-rep(0,bx)

#par(mar=c(5, 8, 4, 6) + 0.1)
plot(log10(Pdiff), axes=F, xlab="",ylab="",type="l",col="black", main=bquote(T[opt]==.(z_vals[species])*degree*C), ylim=c(-6,1))
lines(rep(-4,length(Pdiff)),lty=2, col="black")
axis(4,col="black",lwd=2)
mtext(4,text="Biomass", line=2)
par(new=T)
plot(taudiff, axes=F, xlab="",ylab="",type="l",col="blue", main="")
lines(rep(0,length(taudiff)),lty=2, col="blue")
axis(2,col="blue",lwd=2)
mtext(2,text="Net Dispersal",line=2)
par(new=T)
plot(gmdiff, axes=F, xlab="",ylab="",type="l",col="red", main="")
lines(rep(0,length(gmdiff)),lty=2, col="red")
axis(2,col="red",lwd=2, line=3.5)
mtext(2,text="mu-m",line=5.5)
axis(1,at=1:bx,labels = full_latitude)
mtext("Latitude",side=1,col="black",line=2)

dev.off()


#### Pau ----
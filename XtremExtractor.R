#Library calling
suppressWarnings(suppressMessages(library(ncdf4)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(rnaturalearth)))
suppressWarnings(suppressMessages(library(rnaturalearthdata)))
suppressWarnings(suppressMessages(library(rgeos)))
suppressWarnings(suppressMessages(library(dplyr)))
library(Kendall)
library(biscale)
library(cowplot)
library(ggpubr)
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(raster)
library(modifiedmk)
library(ks)
library(pracma)
library(data.table)
library(RtsEva)
library(xts)
library(regions)

#####Functions#######

mean_daily_value <- function(timeseries) {
  # Convert the 6-hourly time series to daily time series
  daily_timeseries <- apply.daily(timeseries, mean)
  tday=unique(as.Date(as.character(timeseries[,1])))
  return(data.frame(date=tday,Qmd=daily_timeseries))
}

get_POTdata_high <- function(pcts, ms) {
  
  # Calculate the threshold using the specified percentile
  thrsd <- quantile(ms[, 2], pcts / 100)
  
  # Find peaks using the pracma package
  pks_and_locs <- pracma::findpeaks(ms[, 2], minpeakdistance = 7, minpeakheight = thrsd)
  
  # Extract relevant peak information
  pks <- pks_and_locs[, 1]  # Peak values
  locs <- pks_and_locs[, 2] # Indexes of peaks
  st <- pks_and_locs[, 3]   # Start indexes of peaks
  end <- pks_and_locs[, 4]  # End indexes of peaks
  dur <- end-st
  
  f <- ecdf(ms[,2])
  pkq=round(f(pks),5)
  
  # Create the data frame directly with the extracted data
  POTdata_df <- data.frame(
    threshold = thrsd,         # Same threshold value for all rows
    percentile = pcts,         # Same percentile value for all rows
    peaks = pkq,               # Peak values
    stpeaks = st,              # Start index of each peak
    endpeaks = end,            # End index of each peak
    ipeaks = locs,             # Peak locations (indexes)
    time = ms[locs, 1],         # Corresponding time for each peak
    duration = dur
  )
  
  # Return the resulting data frame
  return(POTdata_df)
}

get_POTdata_sliding_low <- function(pcts, ms) {
  
  # Convert the time column (ms[, 1]) to Date objects if they are not already
  ms[, 1] <- as.Date(ms[, 1], origin = "1970-01-01")
  
  # Extract the day of the year (DOY) for each date
  doy <- as.numeric(format(ms[, 1], "%j"))  # %j gives the day of the year (1-366)
  
  # Initialize a vector to store thresholds for each day of the year (1 to 366)
  thresholds <- numeric(366)
  
  # Compute the threshold for each day of the year using a 30-day sliding window (±15 days)
  for (d in 1:366) {
    # Define a window around the current day (±15 days)
    window_days <- c((d - 15):(d + 15)) %% 366  # Wrap around at the end of the year
    window_days[window_days == 0] <- 366  # Correct for modulus result of 0
    
    # Get all data points that fall within this window
    window_data <- ms[doy %in% window_days, 2]
    
    # Calculate the threshold for this day using the specified percentile
    if (length(window_data) > 1) {
      thresholds[d] <- quantile(window_data, pcts / 100)
    } else {
      thresholds[d] <- NA  # If insufficient data, set threshold to NA
    }
  }
  
  # Create a data frame of thresholds for each timestep in ms
  threshold_df <- data.frame(
    time = ms[, 1],  # Time from the input data
    threshold = thresholds[doy],  # Use the thresholds corresponding to the day of the year for each point in time
    measurement = ms[, 2]  # The original measurements
  )
  
  #smoothing of the threshold
  threshold_df$threshold=tsEvaNanRunningMean(threshold_df$threshold,30)
  thresholds= threshold_df$threshold[c(1:366)]
  
  # Initialize a list to store low peak regions
  low_peak_clusters <- list()
  
  # Variables to track cluster start and end
  in_cluster <- FALSE
  cluster_start <- NULL
  cluster_end <- NULL
  
  # Loop over the data to detect low peak clusters
  for (i in seq_len(nrow(ms))) {
    # Get the day of the year for the current data point
    current_doy <- doy[i]
    
    # Get the threshold for this day
    thrsd <- thresholds[current_doy]
    
    # Skip if the threshold is NA (insufficient data)
    if (is.na(thrsd)) next
    
    # Check if the current value is below the dynamic threshold
    if (ms[i, 2] < thrsd) {
      if (!in_cluster) {
        # Start of a new cluster
        in_cluster <- TRUE
        cluster_start <- i
      }
      # Continuously update the cluster end as long as we're below the threshold
      cluster_end <- i
    } else {
      # If we're exiting a low region and were in a cluster, close it
      if (in_cluster) {
        in_cluster <- FALSE
        low_peak_clusters[[length(low_peak_clusters) + 1]] <- c(cluster_start, cluster_end)
        cluster_start <- NULL
        cluster_end <- NULL
      }
    }
  }
  
  # Handle any remaining cluster at the end of the loop
  if (in_cluster) {
    low_peak_clusters[[length(low_peak_clusters) + 1]] <- c(cluster_start, cluster_end)
  }
  
  # Now merge clusters if they are too close (< 30 days apart) and count the number of merged subclusters
  merged_clusters <- list()
  merged_subclusters <- list()
  
  if (length(low_peak_clusters) > 1) {
    current_cluster <- low_peak_clusters[[1]]
    merge_count <- 1  # Start with 1 subcluster
    
    for (j in 2:length(low_peak_clusters)) {
      next_cluster <- low_peak_clusters[[j]]
      
      # Check if the next cluster is within 30 days of the current one
      if (ms[next_cluster[1], 1] - ms[current_cluster[2], 1] <= 30) {
        # Merge the clusters by extending the current cluster's end
        current_cluster[2] <- next_cluster[2]
        merge_count <- merge_count + 1  # Count this as a merged subcluster
      } else {
        # No merge, store the current cluster and move to the next one
        merged_clusters[[length(merged_clusters) + 1]] <- current_cluster
        merged_subclusters[[length(merged_subclusters) + 1]] <- merge_count
        
        # Reset for the next cluster
        current_cluster <- next_cluster
        merge_count <- 1  # Reset merge count for the new cluster
      }
    }
    # Add the last cluster
    merged_clusters[[length(merged_clusters) + 1]] <- current_cluster
    merged_subclusters[[length(merged_subclusters) + 1]] <- merge_count
  } else {
    merged_clusters <- low_peak_clusters
    merged_subclusters <- rep(1, length(low_peak_clusters))  # No merging, each cluster is 1 subcluster
  }
  
  # Prepare the final data frame to store the low peak clusters
  if (length(merged_clusters) > 0) {
    POTdata_df <- data.frame(
      threshold = numeric(),
      percentile = numeric(),
      low_peaks = numeric(),
      stpeaks = numeric(),
      endpeaks = numeric(),
      ipeaks = numeric(),
      time = as.Date(character()),
      cluster_length = numeric(),  # Length of the cluster in days
      merged_subclusters = numeric()  # Number of merged subclusters
    )
    
    # Populate the data frame with the merged low peak clusters
    for (k in seq_along(merged_clusters)) {
      cluster <- merged_clusters[[k]]
      start_idx <- cluster[1]
      end_idx <- cluster[2]
      ipeaks <- which.min(ms[start_idx:end_idx, 2]) + start_idx - 1  # Index of minimum in the cluster
      cluster_length <- as.numeric(ms[end_idx, 1] - ms[start_idx, 1]) + 1  # Length of the cluster
      tw=seq(ms[ipeaks, 1]-15,ms[ipeaks, 1]+15,by="days")
      dotw <- as.numeric(format(tw, "%j")) 
      subset <- ms[which(!is.na(match(doy,dotw))),2]
      f <- ecdf(subset)
      pkf=f(subset)
      pkq=round(f(ms[ipeaks, 2]),5)
      
      POTdata_df <- rbind(POTdata_df, data.frame(
        threshold = thresholds[doy[ipeaks]],
        percentile = pcts,
        low_peaks = pkq,  # The lowest point in the cluster
        stpeaks = ms[start_idx, 1],  # Start date of the cluster
        endpeaks = ms[end_idx, 1],   # End date of the cluster
        ipeaks = ipeaks,             # Index of the lowest point in the cluster
        time = ms[ipeaks, 1],        # Time of the lowest point in the cluster
        cluster_length = cluster_length,  # Length of the cluster in days
        merged_subclusters = merged_subclusters[[k]]  # Number of merged subclusters
      ))
    }
  } else {
    POTdata_df <- data.frame(
      threshold = numeric(),
      percentile = numeric(),
      low_peaks = numeric(),
      stpeaks = numeric(),
      endpeaks = numeric(),
      ipeaks = numeric(),
      time = as.Date(character()),
      cluster_length = numeric(),
      merged_subclusters = numeric()
    )
  }
  
  # Return both the low peak clusters and the threshold data frame
  return(list(low_peaks_df = POTdata_df, threshold_df = threshold_df))
}

outletopen=function(dir,outletname,nrspace=rep(NA,5)){
  ncbassin=paste0(dir,"/",outletname,".nc")
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  if ("Band1"%in% name.vb)namev="Band1"
  name.lon="lon"
  name.lat="lat"
  if (!is.na(nrspace[1])){
    start=as.numeric(nrspace[c(2,4)])
    count=as.numeric(nrspace[c(3,5)])-start+1
  }else{
    londat = ncvar_get(ncb,name.lon) 
    llo=length(londat)
    latdat = ncvar_get(ncb,name.lat)
    lla=length(latdat)
    start=c(1,1)
    count=c(llo,lla)
  }
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll=outll[which(!is.na(outlets)),]
  outlets=outlets[which(!is.na(outlets))]
  outll=data.frame(outlets,outll)
  return (outll)
}

UpAopen=function(dir,outletname,Sloc_final){
  ncbassin=paste0(dir,outletname)
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[2]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncb,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat)
  lla=length(latdat)
  start=c(1,1)
  count=c(llo,lla)
  
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)/1000000
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$upa=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  #outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outll$latlong=paste(round(outll$Var1,4),round(outll$Var2,4),sep=" ")
  outfinal=inner_join(outll, Sloc_final, by="latlong")
  return (outfinal)
}

ReservoirOpen=function(dir,outletname,Sloc_final){
  ncbassin=paste0(dir,outletname)
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncb,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat)
  lla=length(latdat)
  start=c(1,1)
  count=c(llo,lla)
  
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)/1000000
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$upa=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  #outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outll$latlong=paste(round(outll$Var1,4),round(outll$Var2,4),sep=" ")
  outfinal=inner_join(outll, Sloc_final, by="latlong")
  return (outfinal)
}


###### Data loading @@@@@@@@@@@
#Set data directory
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")

#load outf

outf=c()
for( Nsq in 1:88){
  print(Nsq)
  rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
  rspace=rspace[,-1]
  nrspace=rspace[Nsq,]
  #outletname="outletsv8_hybas07_01min"
  #outletname="outlets_hybas09_01min"
  outletname="efas_rnet_100km_01min"
  
  outhybas=outletopen(hydroDir,outletname,nrspace)
  Idstart=as.numeric(Nsq)*10000
  Idstart2=as.numeric(Nsq)*100000
  if (length(outhybas$outlets)>0){
    outhybas$outlets=seq((Idstart+1),(Idstart+length(outhybas$outlets)))
    outhybas$outl2=seq((Idstart2+1),(Idstart2+length(outhybas$outlets)))
    outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")
    #outcut=which(!is.na(match(outhybas$outlets,parlist$catchment)))
    # zebi=seq(parlist$catchment[1],parlist$catchment[length(parlist$catchment)])
    # outcut=which(!is.na(match(outhybas$outlets,zebi)))
    outhloc=outhybas
    outf=rbind(outf,outhloc)
  }
}


### HydroRegions ----

# GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
# GHR=as.data.frame(GridHR,xy=T)
# GHR=GHR[which(!is.na(GHR[,3])),]
# GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ") 
# GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
# GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted.shp")
# HydroRsf=fortify(GHshpp) 


### NUTS3 ----
# NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_Extended_domain.shp"))
# NUTS3$N3ID=c(1:length(NUTS3$NUTS_ID))
# N2ID=unique(NUTS3$NUTS2_ID)
# N2IDn=c(1:length(N2ID))
# mati=match(NUTS3$NUTS2_ID,N2ID)
# NUTS3$N2ID=N2IDn[mati]
#st_write(NUTS3, paste0(hydroDir,"/Countries/NUTS3/NUTS3_modified.shp"), driver = "ESRI Shapefile")

NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_modified.shp"))
GridNUTS3=raster( paste0(hydroDir,"/Countries/NUTS3/NUTS3_Raster3ID.tif"))
GN3=as.data.frame(GridNUTS3,xy=T)
GN3=GN3[which(!is.na(GN3[,3])),]
GN3$llcoord=paste(round(GN3$x,4),round(GN3$y,4),sep=" ")
GN3_riv=right_join(GN3,outf,by= c("llcoord"="latlong"))

GridNUTS2=raster( paste0(hydroDir,"/Countries/NUTS3/NUTS3_Raster2ID.tif"))
GN2=as.data.frame(GridNUTS2,xy=T)
GN2=GN2[which(!is.na(GN2[,3])),]
GN2$llcoord=paste(round(GN2$x,4),round(GN2$y,4),sep=" ")
GN2_riv=right_join(GN2,outf,by= c("llcoord"="latlong"))

GNF=right_join(GN3,GN2_riv,by="llcoord")

GNUTS3sf=fortify(NUTS3)

GNFx=GNF[which(is.na(GNF$NUTS3_Raster3ID)),]

#Load old and new NUTS3 region IDs
NUTS3_2010 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/Regions_v2010_simplified.shp"))

NUTS3_2021 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/Regions_v2021_simplified.shp"))


cores=nuts_changes[which(nuts_changes$typology=="nuts_level_3"),]

#load UpArea
#load upstream area
main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
outletname="/GIS/upArea_European_01min.nc"
dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(valid_path,outletname,outf)
head(UpArea)

###Plot parameters###
palet2=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = TRUE, fixup = TRUE))
outletname="efas_rnet_100km_01min"
outll=outletopen(hydroDir,outletname)
cord.dec=outll[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
Impdates=seq(1950,2020,by=10)
valuenames=paste0("Y",Impdates)
catmap=cst7
basemap=w2



###Generate daily Q#####
pcts=95
yseq=seq(1980,2020)
totime=c()
NR_Qd_total=c()
for (year in yseq){
  print(year)
  file_name <- paste0("/NUTSQ/NUTS3x_", year, ".csv")
  NUTsQ<-read.csv(file=paste0(hydroDir,file_name))
  
  #extract the list of NUTS regions
  
  Nreg=unique(NUTsQ$V1)
  
  #match these ID with real NUTS ID
  mn=match(Nreg,NUTS3$N3ID)
  Nutlist=NUTS3[mn,]
  
  NR1=NUTsQ[which(NUTsQ$V1==Nreg[1]),]
  NRtime=as.POSIXct(NR1$V3*24*3600-3600,origin="1979-01-01 00:00:00")
  NR_Qd=c()
  for (n in 1:length(Nreg)){
    print(n)
    NR1=NUTsQ[which(NUTsQ$V1==Nreg[n]),]
    #daily aggregation
    NR_short=data.frame(time=NRtime, Q=NR1$V2)
    NR_day=mean_daily_value(NR_short)
    NR_Qd=cbind(NR_Qd,NR_day$Q)
    #plot(NR_day,type="o")
  }
  timevec=NR_day$date
  totime=c(totime,timevec)
  NR_Qd_total=rbind(NR_Qd_total,NR_Qd)
}


###load daily Q####
NR_Qd_total<-read.csv(file=paste0(hydroDir,"Daily_Q_1980-2020.csv"))
NR_Qd_total=NR_Qd_total[,-1]
NRtime=seq(1:length(NR_Qd_total$HR064))
NRtime=as.Date(NRtime,origin="1979-12-31")
NRtime=as.POSIXct(NRtime)



####load Dominik's event set to do a first matching####
Hanze_flood=read.csv(file="D:/tilloal/Documents/01_Projects/RiskDynamics/Data/HANZE_events.csv")
#I need to create a vector of affected regions for each event
list_of_word_vectors <- lapply(Hanze_flood$Regions.affected..v2021., function(x) unlist(strsplit(x, ";")))
vectors_2021=unlist(list_of_word_vectors)
vectors_2021=unique(vectors_2021)
list_of_word_vectors2 <- lapply(Hanze_flood$Regions.affected..v2010., function(x) unlist(strsplit(x, ";")))

vectors_2010=unlist(list_of_word_vectors2)
vectors_2010=unique(vectors_2010)
N1021=cbind(vectors_2010,vectors_2021)
#extract one location

NutVector=NUTS3$NUTS_ID



df_POTh_large=c()
df_POTl_large=c()
Hanze_flood_BN=c()
High_events_df=c()
Low_events_df=c()
for (id in 1:length(NutVector)){
  print(id)
  NUTl=NutVector[id]
  myloc=NR_Qd_total[,id]
  if (length(which(!is.na(myloc)))>0){
    ms=data.frame(time=NRtime,data=myloc)
    
    POTdata_high <- get_POTdata_high(98,ms)
    POTdata_high$enddate=as.POSIXct(POTdata_high$endpeaks*24*3600-3600,origin="1980-01-01 00:00:00")
    POTdata_low <- get_POTdata_sliding_low(05, ms)
    
    #plot(POTdata_low$threshold_df$threshold[1:366], type="l",xlab="day of the year", ylab="Q", ylim=c(0,0.5), lwd=2)
    #saving raw events detected
    Low_events=POTdata_low$low_peaks_df
    Low_events$NUT=NUTl
    High_events=POTdata_high
    High_events$NUT=NUTl
    
    Low_events_df=rbind(Low_events_df,Low_events)
    High_events_df=rbind(High_events_df,High_events)
    
    
    matches <- lapply(list_of_word_vectors, function(vector) NUTl %in% vector)
    tm=which(matches==T)
    Hanze_Nut=Hanze_flood[tm,]
    Hanze_Nut=Hanze_Nut[which(Hanze_Nut$Year>1980),]
    if(length(Hanze_Nut$Country.code)>0){
      Hanze_bNut=data.frame(NUTl,Hanze_Nut)
      Hanze_flood_BN=rbind(Hanze_flood_BN,Hanze_bNut)
    }
    #now check overlap with mt detected anomalies
    df_POTh_keep=c()
    df_POTl_keep=c()
    for (e in 1:length(Hanze_Nut$ï..ID)){
      d1=(as.Date(POTdata_high$enddate)-as.Date(Hanze_Nut$End.date[e]))
      keep=which(d1<7 & d1>-14)
      if (length(keep)>0){
        POTh_keep=POTdata_high[keep,]
        POTh_keep$dtime=d1[keep]
        POTh_keep$NUTID=NUTl
        POTh_keep$eventID=Hanze_Nut$ï..ID[e]
        df_POTh_keep=rbind(df_POTh_keep,POTh_keep)
      } else{
        
      }
      
      d2=(as.Date(POTdata_low$low_peaks_df$endpeaks)-as.Date(Hanze_Nut$End.date[e]))
      keep=which(d2<0 & d2>-14)
      if (length(keep)>0){
        POTl_keep=POTdata_low$low_peaks_df[keep,]
        POTl_keep$dtime=d2[keep]
        POTl_keep$NUTID=NUTl
        POTl_keep$eventID=Hanze_Nut$ï..ID[e]
        df_POTl_keep=rbind(df_POTl_keep,POTl_keep)
      }
    }
    
    df_POTh_keep=as.data.frame(df_POTh_keep)
    df_POTl_keep=as.data.frame(df_POTl_keep)
    
    df_POTh_large=rbind(df_POTh_large,df_POTh_keep)
    df_POTl_large=rbind(df_POTl_large,df_POTl_keep)
  }
  
}

#save outputs from here

save(High_events_df,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/events_highflows_1980-2020_v3.Rdata"))
#save(Low_events_df,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/events_lowflows_1980-2020_v2.Rdata"))

save(df_POTh_large,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/HighflowsXHanze_1980-2020_v3.Rdata"))
#save(df_POTl_large,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/LowflowsXHanze_1980-2020_v2.Rdata"))

#new loop for floods with varying threshold


load(file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/events_highflows_1980-2020.Rdata"))
load(file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/events_lowflows_1980-2020.Rdata"))
load(file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/HighflowsXHanze_1980-2020.Rdata"))
load(file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/LowflowsXHanze_1980-2020.Rdata"))


#save events in csv


write.csv(High_events_df,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Hflow_Q_1980-2020.csv"))
write.csv(Low_events_df,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Lflow_Q_1980-2020.csv"))




#########Aggregate flood and drought events in time##############

#This is in the objective of creating a dataet similar to HANZE

#loop over time ?
# then I have to segregate by distance
length(High_events_df$NUT)/ length(NUTS3$NUTS_ID)/40

#for each day in the daylist, where is there also another event

undate=unique(High_events_df$time)
undate=as.Date(undate)
undate=undate[order(undate)]
zdate=seq.Date(undate[1],undate[length(undate)],by="days")

neighbors_list <- st_touches(GNUTS3sf)
GNUTS3sf$Neigh <- sapply(seq_len(length(neighbors_list)), function(i) {
  neighbor_ids <- GNUTS3sf$NUTS_ID[neighbors_list[[i]]] # Retrieve the NUTS_ID of each neighbor
  paste(neighbor_ids, collapse = " - ")          # Combine the NUTS_IDs into a single string separated by " - "
})                                               # Create a new column 'Neigh' to store the NUTS_ID of neighbors



# Load data.table for faster data manipulation
library(data.table)

# Convert High_events_df to data.table for faster access
High_events_dt <- as.data.table(High_events_df)
setkey(High_events_dt, NUT)

# Precompute neighbors for each NUTS region
neighbor_list <- lapply(GNUTS3sf$NUTS_ID, function(nut_id) {
  nei <- GNUTS3sf$Neigh[match(nut_id, GNUTS3sf$NUTS_ID)]
  unlist(strsplit(nei, " - "))
})
names(neighbor_list) <- GNUTS3sf$NUTS_ID

# Initialize idt_ev as NA
High_events_dt[, idt_ev := NA_character_]

# Loop over each date in undate
for (d in 1:length(undate)) {
  # Define the date and calculate time differences
  print(d)
  day <- undate[d]
  ladif <- abs(as.numeric(difftime(as.Date(High_events_dt$time), day, units = "days")))
  
  # Select rows within the 3-day temporal window
  rowsel <- which(ladif < 3)
  
  trouduc=High_events_dt[rowsel,]
  # For rows within the temporal window, get unique regions
  Nutin <- unique(High_events_dt[rowsel,NUT])
  
  # Track clusters for merged events
  cluster_id <- 1
  for (n in Nutin) {
    # Get neighbors for the current NUTS region
    sel <- neighbor_list[[n]]
    # Find matching rows within the selected regions and their neighbors
    ev_nut <- c(n, intersect(sel, Nutin))
    row_ev <- rowsel[High_events_dt[rowsel, NUT] %in% ev_nut] 
    # Check if any row in the event group already has an idt_ev
    #existing_id <- unique(na.omit(High_events_dt$idt_ev[row_ev]))
    existing_id <- unique(na.omit(High_events_dt[row_ev, idt_ev]))
    #print(existing_id)
    if (length(existing_id) > 0) {
      # Use the existing id if found
      High_events_dt[row_ev, idt_ev := existing_id[1]]
    } else {
      # Assign a new id for this cluster if no existing id is found
      new_id <- paste(d, cluster_id, sep = "-")
      High_events_dt[row_ev, idt_ev := new_id]
      cluster_id <- cluster_id + 1
    }
  }
}

#####################################################

#save(High_events_dt,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/Modelled_events_1980-2020.Rdata"))

load(file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/Modelled_events_1980-2020.Rdata"))
ev_unik=(unique(High_events_dt$idt_ev))
ev_uid=c(1:length(ev_unik))

#aggregatiion

event_modelled=aggregate(list(events=High_events_dt$duration),
                   by = list(HYBAS_ID=High_events_dt$idt_ev),
                   FUN = function(x) c(mean=mean(x,na.rm=T),len=length(x)))
event_modelled <- do.call(data.frame, event_modelled)

event_modelled$id=ev_uid

is_integer_value <- function(val) {
  val == as.integer(val)
}

lem=length(event_modelled$id)
#I would need another loop to identify beginning and end of each event
st_vec=vector("numeric",length=lem)
en_vec=vector("numeric",length=lem)
Anut_vec=vector("character",length=lem)
for (eid in 1:lem)
{
  if (is_integer_value(eid/1000)) print(eid)
  evd=event_modelled$HYBAS_ID[eid]
  
  ziz=which(High_events_dt$idt_ev==evd)
  myev=High_events_dt[ziz,]
  stdate=min(myev$stpeaks)
  edate=max(myev$endpeaks)
  
  involved_NUTS=myev$NUT
  Anuts <- paste(involved_NUTS, collapse = "-")
  Anut_vec[eid]= Anuts
  st_vec[eid]=stdate
  en_vec[eid]=edate
}

# Convert back to data.frame if needed
event_modelled <- data.frame(event_modelled,st_vec,en_vec,Anut_vec)

save(event_modelled,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/Modelled_floodevents_grouped_1980-2020.Rdata"))


##################################################################################################

#Load Dominik's corrected impact database
#load Dominik's event set to do a first matching
Hanze_events_1950=read.csv(file="D:/tilloal/Documents/01_Projects/RiskDynamics/Data/Attribution_losses_under_1950_exposure_protection_and_vulnerability.csv")
Hanze_events_1950_N2021=left_join(Hanze_events_1950,cores,by=c("NUTS3"="code_2010"))
Hanze_events_1950_N2021$code_2021[which(is.na(Hanze_events_1950_N2021$code_2021))]=Hanze_events_1950_N2021$NUTS3[which(is.na(Hanze_events_1950_N2021$code_2021))]
mh=match(Hanze_events_1950$Index,Hanze_events_1950_N2021$Index)
Hanze_events_1950$NUTS2021=Hanze_events_1950_N2021$code_2021[mh]
#Matching Domi's events with detected flood events (same events)

NutVector=NUTS3$NUTS_ID
Hanze_flood_BN=c()
High_eventsxHanze=c()
Low_eventsxHanze=c()
NUTno=c()
NUTnow=c()
for (id in 1:length(NutVector)){
  print(id)
  NUTl=NutVector[id]

  Low_events=Low_events_df[which(Low_events_df$NUT==NUTl),]
  High_events=High_events_df[which(High_events_df$NUT==NUTl),]
  
  Hanze_Nut=Hanze_events_1950[which(Hanze_events_1950$NUTS2021==NUTl),]
  Hanze_Nut=Hanze_Nut[which(Hanze_Nut$Year>=1980),]
  if(length(Hanze_Nut$Year)==0){
    NUTno=c(NUTno,NUTl)
  }
  matches <- lapply(list_of_word_vectors, function(vector) NUTl %in% vector)
  tm=which(matches==T)
  Hanze_Nut2=Hanze_flood[tm,]
  Hanze_Nut2=Hanze_Nut2[which(Hanze_Nut2$Year>1980),]
  
  if(length(Hanze_Nut2$Year)==0){
    NUTnow=c(NUTnow,NUTl)
  }
  if(length(Hanze_Nut2$Year)>0){
    Hanze_Nut2$NUTl=NUTl
    Hanze_flood_BN=rbind(Hanze_flood_BN,Hanze_Nut2)
  }
  #now check overlap with mt detected anomalies
  df_POTh_keep=c()
  df_POTl_keep=c()
  for (e in 1:length(Hanze_Nut2$ï..ID)){
    d1=(as.Date(High_events$enddate)-as.Date(Hanze_Nut$End_date[e]))
    keep=which(d1<7 & d1>(-90))
    if (length(keep)>0){
      POTh_keep=High_events[keep,]
      POTh_keep$dtime=d1[keep]
      POTh_keep$NUTID=NUTl
      POTh_keep$eventID=Hanze_Nut2$ï..ID[e]
      High_eventsxHanze=rbind(High_eventsxHanze,POTh_keep)
    } else{
      
    }
    
    
    d2=(as.Date(Low_events$endpeaks)-as.Date(Hanze_Nut$End_date[e]))
    keep=which(d2<0 & d2>-90)
    if (length(keep)>0){
      POTl_keep=Low_events[keep,]
      POTl_keep$dtime=d2[keep]
      POTl_keep$NUTID=NUTl
      POTl_keep$eventID=Hanze_Nut2$ï..ID[e]
      Low_eventsxHanze=rbind(Low_eventsxHanze,POTl_keep)
    }
  }
  
  High_eventsxHanze=as.data.frame(High_eventsxHanze)
  Low_eventsxHanze=as.data.frame(Low_eventsxHanze)
}

length(unique(Hanze_flood_BN$ï..ID))
Hanze_flood_1980=Hanze_flood[which(Hanze_flood$Year>1980),]

ptt=match(Hanze_flood_1980$ï..ID,unique(Hanze_flood_BN$ï..ID))
arg=which(is.na(ptt))
Hanze_flood_1980[arg,]
#Matching Domi's events with previous flood events (antecedent condition is wet)
idfd=which(!is.na(match(Hanze_events_1950$NUTS2021,NutVector)))
Hanze_events_1950_fd=Hanze_events_1950[idfd,]
Hanze_events_1980_fd=Hanze_events_1950[which(Hanze_events_1950_fd$Year>1980),]

High_eventsxHanze_m=High_eventsxHanze[which(abs(High_eventsxHanze$dtime)<14),]

length(High_eventsxHanze_m$threshold)/length(Hanze_flood_BN$Country.code)

#Looking for discrapancies in modelled damages

#load csv file of modelled damages

modelled_events=read.csv(file="D:/tilloal/Documents/01_Projects/RiskDynamics/Data/modelled/HANZE_potential_flood_catalogue_all.csv")
modelled_events$HANZE.ID=as.numeric(modelled_events$HANZE.ID)
length(which(!is.na(modelled_events$HANZE.ID)))
Hanze_obsmod=left_join(Hanze_flood,modelled_events,by=c("ï..ID"="HANZE.ID"))
#Matching Domi's events with previous drought events (antecedent condition is dry)


###################################################################################
#Cascading scenarios...


mynut="FRJ12"
Reg_HANZE=Hanze_flood_BN[ which(Hanze_flood_BN$NUTl=="FRJ12"),]
Reg_Flood=High_events_df[which(High_events_df$NUT=="FRJ12"),]
Reg_Drought=Low_events_df[which(Low_events_df$NUT=="FRJ12"),]

Reg_HANZE$y_fixed=1
Reg_Flood$y_fixed=2
Reg_Drought$y_fixed=3
Reg_HANZE$time=as.Date(Reg_HANZE$Start.date)
palet=c(hcl.colors(9, palette = "Reds", alpha = NULL, rev = T, fixup = TRUE))
ggplot() +
  geom_point(data=Reg_HANZE, aes(x = time, y = y_fixed, size = Fatalities,fill=Losses..2020.euro.),shape=21, col="black") +
  geom_point(data=Reg_Flood, aes(x = as.Date(time), y = y_fixed),shape=24,fill="royalblue",col="black",size=5) +
  geom_point(data=Reg_Drought, aes(x = time, y = y_fixed),shape=25,fill="orange",size=5) +
  scale_y_continuous(limits=c(0,4),breaks=c(1,2,3),labels=c("HANZE floods","High flow anomalies","low flow anomalies"))+
  scale_fill_gradientn(
    colors=palet, name="Economic losses\n(Million euros)",breaks=c(1e7,1e8,1e9),trans="log",labels=c(10,100,1000),limits=c(1e7,2e9))+ 
  # Plot points with varying size
  scale_size_continuous(range = c(1, 15)) +  # Set point size range
  labs(x = "Time", y = "Events") +  # Axis labels
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major.y  = element_line(colour = "black"),
        panel.grid.major.x  = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle(paste0("Timeline for NUTS3 region ",mynut))


ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/timeline_plot.jpg"), width=40, height=16, units=c("cm"),dpi=300) 




##############################################

#Create a fast loop exploiting data previously generated for different scenarii

#load events from Michele

hw_events=read.csv(file="D:/tilloal/Documents/01_Projects/RiskDynamics/Data/multihz_data/hw_nuts3_v1.csv")

shortHW=c(7,-7)
mediumHW=c(7,-14)


shortFlood=c(7,-14)
mediumFlood=c(7,-30)
longFlood=c(7,-60)

shortDrought=c(0,-14)
mediumDrought=c(0,-30)
longDrought=c(0,-90)


NutVector=NUTS3$NUTS_ID

wflood=shortFlood
wdrought=shortDrought
wheat=shortHW
df_POTx_large=c()
df_POTh_large=c()
df_POTl_large=c()
Hanze_flood_BN=c()
for (id in 1:length(NutVector)){
  print(id)
  NUTl=NutVector[id]
  myloc=NR_Qd_total[,id]
  if (length(which(!is.na(myloc)))>0){
    ms=data.frame(time=NRtime,data=myloc)
    Low_events=Low_events_df[which(Low_events_df$NUT==NUTl),]
    High_events=High_events_df[which(High_events_df$NUT==NUTl),]
    HW_events=hw_events[which(hw_events$nuts_id==NUTl),]
    matches <- lapply(list_of_word_vectors, function(vector) NUTl %in% vector)
    tm=which(matches==T)
    Hanze_Nut=Hanze_flood[tm,]
    Hanze_Nut=Hanze_Nut[which(Hanze_Nut$Year>1980),]
    if(length(Hanze_Nut$Country.code)>0){
      Hanze_bNut=data.frame(NUTl,Hanze_Nut)
      Hanze_flood_BN=rbind(Hanze_flood_BN,Hanze_bNut)
    }
    #now check overlap with mt detected anomalies
    df_POTh_keep=c()
    df_POTl_keep=c()
    df_POTx_keep=c()
    for (e in 1:length(Hanze_Nut$ï..ID)){
      d1=(as.Date(High_events$enddate)-as.Date(Hanze_Nut$End.date[e]))
      keep=which(d1<wflood[1] & d1>wflood[2])
      if (length(keep)>0){
        POTh_keep=High_events[keep,]
        POTh_keep$dtime=d1[keep]
        POTh_keep$NUTID=NUTl
        POTh_keep$eventID=Hanze_Nut$ï..ID[e]
        df_POTh_keep=rbind(df_POTh_keep,POTh_keep)
      } else{
        
      }
    
    d2=(as.Date(Low_events$endpeaks)-as.Date(Hanze_Nut$End.date[e]))
    keep=which(d2<wdrought[1] & d2>wdrought[2])
    if (length(keep)>0){
      POTl_keep=Low_events[keep,]
      POTl_keep$dtime=d2[keep]
      POTl_keep$NUTID=NUTl
      POTl_keep$eventID=Hanze_Nut$ï..ID[e]
      df_POTl_keep=rbind(df_POTl_keep,POTl_keep)
    }
    
    d3=(as.Date(HW_events$end)-as.Date(Hanze_Nut$End.date[e]))
    keep=which(d3<wheat[1] & d3>wheat[2])
    if (length(keep)>0){
      POTx_keep=HW_events[keep,]
      POTx_keep$dtime=d2[keep]
      POTx_keep$NUTID=NUTl
      POTx_keep$eventID=Hanze_Nut$ï..ID[e]
      df_POTx_keep=rbind(df_POTx_keep,POTx_keep)
    }
  }
  
  df_POTh_keep=as.data.frame(df_POTh_keep)
  df_POTl_keep=as.data.frame(df_POTl_keep)
  df_POTx_keep=as.data.frame(df_POTx_keep)
  
  df_POTh_large=rbind(df_POTh_large,df_POTh_keep)
  df_POTl_large=rbind(df_POTl_large,df_POTl_keep)
  df_POTx_large=rbind(df_POTx_large,df_POTx_keep)
  
  }
}

#ok good stuff
# I can also decide to limit lag to one month
Hanze_flood_BN$NEID=paste0(Hanze_flood_BN$ï..ID,Hanze_flood_BN$NUTl)
Hanze_flood_BN$NEID
df_POTh_large$NEID=paste0(df_POTh_large$eventID,df_POTh_large$NUTID)

Hanze_flood_1980=Hanze_flood[which(Hanze_flood$Year>1980),]
#Now matching with all events to get basic metrics

Events_prefflood=inner_join(df_POTh_large,Hanze_flood_1980, by=c("eventID"="ï..ID"))

#Events_prefflood=full_join(df_POTh_large,Hanze_flood_BN, by=c("NEID"))

Events_precoflood=inner_join(df_POTl_large,Hanze_flood_1980, by=c("eventID"="ï..ID"))

Events_prehflood=inner_join(df_POTx_large,Hanze_flood_1980, by=c("eventID"="ï..ID"))

#verification of missing events
v1=(unique(Hanze_flood_1980$ï..ID))
v2=(unique(Events_prefflood$ï..ID))
vcrap=v1[which(is.na(match(v1,v2)))]
oyoy=Hanze_flood_1980[match(vcrap,Hanze_flood_1980$ï..ID),]



#aggregate by event ID
Ev_flmatch=aggregate(list(time=Events_prefflood$time),
                    by = list(NID=Events_prefflood$NUTID,EID=Events_prefflood$eventID),
                    FUN = function(x) c(l=length(x)))
names(Ev_flmatch)[3]="N_PreEv"
Ev_flmatch$NEID=paste0(Ev_flmatch$EID,Ev_flmatch$NID)



Ev_fli=aggregate(list(time=Ev_flmatch$N_PreEv),
                 by = list(EID=Ev_flmatch$EID),
                 FUN = function(x) c(l=length(x)))

Ev_flF=aggregate(list(time=Ev_floods$N_PreEv),
                 by = list(EID=Ev_floods$EID),
                 FUN = function(x) c(l=length(x)))

length(Ev_fli$EID[which(Ev_fli$time>1)])

Hanze_flood_BN$NEID=paste0(Hanze_flood_BN$ï..ID,Hanze_flood_BN$NUTl)

Ev_flood2=full_join(Ev_flmatch,Hanze_flood_BN,by="NEID")

Ev_flood2$N_PreEv[which(is.na(Ev_flood2$N_PreEv))]=0
Ev_flood2$NID[which(is.na(Ev_flood2$NID))]=Ev_flood2$NUTl[which(is.na(Ev_flood2$EID))]
Ev_flood2$EID[which(is.na(Ev_flood2$EID))]=Ev_flood2$ï..ID[which(is.na(Ev_flood2$EID))]

#Identify events for which the previous event is the event itself

Ev_flood2$PC=Ev_flood2$N_PreEv
Ev_flood2$PC[which(Ev_flood2$PC<2)]=0
Ev_flood2$PC[which(Ev_flood2$PC>1)]=1


# Important metrics
Ev_flF=aggregate(list(time=Ev_flood2$PC),
                 by = list(EID=Ev_flood2$EID),
                 FUN = function(x) c(l=mean(x)))

#percentage of CE
a=length(Ev_flF$time[which(Ev_flF$time==0)])
b=length(Ev_flF$time[which(Ev_flF$time>0)])
b/(a+b)*100

Ev_flX=aggregate(list(time=Ev_flood2$Losses..2020.euro.),
                 by = list(EID=Ev_flood2$EID),
                 FUN = function(x) c(l=mean(x)))

#Loss ratio
a=sum(Ev_flX$time[which(Ev_flF$time==0)],na.rm=T)
b=sum(Ev_flX$time[which(Ev_flF$time>0)],na.rm=T)
b/(a)


#same for drought



#aggregate by event ID
Ev_drmatch=aggregate(list(time=Events_precoflood$time),
                     by = list(NID=Events_precoflood$NUTID,EID=Events_precoflood$eventID),
                     FUN = function(x) c(l=length(x)))
names(Ev_drmatch)[3]="N_PreEv"
Ev_drmatch$NEID=paste0(Ev_drmatch$EID,Ev_drmatch$NID)



Ev_dri=aggregate(list(time=Ev_drmatch$N_PreEv),
                 by = list(EID=Ev_drmatch$EID),
                 FUN = function(x) c(l=length(x)))


length(Ev_dri$EID[which(Ev_dri$time>1)])

Hanze_flood_BN$NEID=paste0(Hanze_flood_BN$ï..ID,Hanze_flood_BN$NUTl)

Ev_drought2=full_join(Ev_drmatch,Hanze_flood_BN,by="NEID")

Ev_drought2$N_PreEv[which(is.na(Ev_drought2$N_PreEv))]=0
Ev_drought2$NID[which(is.na(Ev_drought2$NID))]=Ev_drought2$NUTl[which(is.na(Ev_drought2$EID))]
Ev_drought2$EID[which(is.na(Ev_drought2$EID))]=Ev_drought2$ï..ID[which(is.na(Ev_drought2$EID))]

#Identify events for which the previous event is the event itself

Ev_drought2$PC=Ev_drought2$N_PreEv
Ev_drought2$PC[which(Ev_drought2$PC<1)]=0
Ev_drought2$PC[which(Ev_drought2$PC>=1)]=1


# Important metrics
Ev_flF=aggregate(list(time=Ev_drought2$PC),
                 by = list(EID=Ev_drought2$EID),
                 FUN = function(x) c(l=mean(x)))

#percentage of CE
a=length(Ev_flF$time[which(Ev_flF$time==0)])
b=length(Ev_flF$time[which(Ev_flF$time>0)])
b/(a+b)*100

Ev_flX=aggregate(list(time=Ev_drought2$Losses..2020.euro.),
                 by = list(EID=Ev_drought2$EID),
                 FUN = function(x) c(l=mean(x)))

#Loss ratio
a=sum(Ev_flX$time[which(Ev_flF$time==0)],na.rm=T)
b=sum(Ev_flX$time[which(Ev_flF$time>0)],na.rm=T)
b/(a)


#same for HW



#aggregate by event ID
Ev_Xmatch=aggregate(list(time=Events_prehflood$begin),
                     by = list(NID=Events_prehflood$NUTID,EID=Events_prehflood$eventID),
                     FUN = function(x) c(l=length(x)))
names(Ev_Xmatch)[3]="N_PreEv"
Ev_Xmatch$NEID=paste0(Ev_Xmatch$EID,Ev_Xmatch$NID)



Ev_Xi=aggregate(list(time=Ev_Xmatch$N_PreEv),
                 by = list(EID=Ev_Xmatch$EID),
                 FUN = function(x) c(l=length(x)))


length(Ev_Xi$EID[which(Ev_Xi$time>1)])

Hanze_flood_BN$NEID=paste0(Hanze_flood_BN$ï..ID,Hanze_flood_BN$NUTl)

Ev_heat2=full_join(Ev_Xmatch,Hanze_flood_BN,by="NEID")

Ev_heat2$N_PreEv[which(is.na(Ev_heat2$N_PreEv))]=0
Ev_heat2$NID[which(is.na(Ev_heat2$NID))]=Ev_heat2$NUTl[which(is.na(Ev_heat2$EID))]
Ev_heat2$EID[which(is.na(Ev_heat2$EID))]=Ev_heat2$ï..ID[which(is.na(Ev_heat2$EID))]

#Identify events for which the previous event is the event itself

Ev_heat2$PC=Ev_heat2$N_PreEv
Ev_heat2$PC[which(Ev_heat2$PC<1)]=0
Ev_heat2$PC[which(Ev_heat2$PC>=1)]=1


# Important metrics
Ev_flF=aggregate(list(time=Ev_heat2$PC),
                 by = list(EID=Ev_heat2$EID),
                 FUN = function(x) c(l=mean(x)))

#percentage of CE
a=length(Ev_flF$time[which(Ev_flF$time==0)])
b=length(Ev_flF$time[which(Ev_flF$time>0)])
b/(a+b)*100

Ev_flX=aggregate(list(time=Ev_heat2$Losses..2020.euro.),
                 by = list(EID=Ev_heat2$EID),
                 FUN = function(x) c(l=mean(x)))

#Loss ratio
a=sum(Ev_flX$time[which(Ev_flF$time==0)],na.rm=T)
b=sum(Ev_flX$time[which(Ev_flF$time>0)],na.rm=T)
b/(a)



##############A deep cleaning is required here




Ev_single=Ev_flood2[which(Ev_flood2$PC==0),]
length(unique(Ev_single$EID))
Ev_compound=Ev_flood2[which(Ev_flood2$PC==1),]
length(unique(Ev_compound$EID))

length(Ev_flood2$NID.x)
Regio_mplot=aggregate(list(time=Ev_flood2$Losses..2020.euro.),
                        by = list(NID=Ev_flood2$NID.x,Preflood=Ev_flood2$PC),
                        FUN = function(x) c(l=length(x),mean=mean(x,na.rm=T),median=median(x,na.rm=T)))
Regio_mplot <- do.call(data.frame, Regio_mplot)

Regio_mplot0=Regio_mplot[which(Regio_mplot$Preflood==0),]
Regio_mplot1=Regio_mplot[which(Regio_mplot$Preflood==1),]

length(Regio_mplot0$NID)
length(Regio_mplot1$NID)
Regio_mplot0=inner_join(Regio_mplot1,Regio_mplot0,by="NID")
Regio_mplot0$Ratio=Regio_mplot0$time.mean.x/Regio_mplot0$time.mean.y
#NUTS3 maps
Regio_mplot=full_join(Regio_mplot0,NUTS3,by=c("NID"="NUTS_ID"))

#Nplot <- st_transform(Regio_mplot, crs = 3035)


nuts1 <- st_read("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/NUTS_RG_20M_2021_3035.shp")


Regio_mplot=st_as_sf(Regio_mplot)
nuts1=nuts1[which(nuts1$LEVL_CODE==1),]
n3to1 <- st_join( nuts1, Regio_mplot,join = st_intersects)


#aggregate byy nuts1


Regio_n1plot_cp=aggregate(list(val=n3to1$time.mean.x),
                          by = list(NID=n3to1$NUTS_ID),
                          FUN = function(x) c(l=length(x),sum=sum(x,na.rm=T),mean=mean(x,na.rm=T),median=median(x,na.rm=T)))
Regio_n1plot_cp <- do.call(data.frame, Regio_n1plot_cp)


Regio_n1plot_sg=aggregate(list(val=n3to1$time.mean.y),
                          by = list(NID=n3to1$NUTS_ID),
                          FUN = function(x) c(l=length(x),sum=sum(x,na.rm=T),mean=mean(x,na.rm=T),median=median(x,na.rm=T)))
Regio_n1plot_sg <- do.call(data.frame, Regio_n1plot_sg)

Regio_n1plot_cp$ratio=(Regio_n1plot_cp$val.sum+1)/(Regio_n1plot_sg$val.sum+1)
Regio_n1plot_cp$ratio[which(Regio_n1plot_cp$val.sum==0 | Regio_n1plot_sg$val.sum==0)]=NA


Regio_n1plot=full_join(Regio_n1plot_cp,nuts1,by=c("NID"="NUTS_ID"))


palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = T, fixup = TRUE))

pl1=ggplot(basemap) +
  geom_sf(fill="gray95")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=Regio_n1plot,aes(fill=ratio,geometry=geometry),color="gray12",alpha=1,size=0.25,stroke=0,shape=15)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=paletf,limits=c(0.1,10),trans="log", breaks=c(1e-2,0.1,1,10,100),
    oob = scales::squish,na.value="gray30", name="Loss ratio") +
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("Ratio of losses from preconditionned flood events over single flood events (high flows)")
pl1


ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/LossRatioMap_hflow.jpg"),pl1, width=30, height=25, units=c("cm"),dpi=800) 







Hanze_floodEV=aggregate(list(time=Hanze_flood_BN$Year),
                        by = list(NID=Hanze_flood_BN$NUTl),
                        FUN = function(x) c(l=length(x)))



Ev_flood3=aggregate(list(Ev_flood2$precon),
                    by = list(EID=Ev_flood2$EID.x),
                    FUN = function(x) c(l=length(x),s=mean(x)))
Ev_flood3 <- do.call(data.frame, Ev_flood3)
names(Ev_flood3)[c(2,3)]=c("Nregion","NpreEv")

length(Ev_flood3$EID)/length(Hanze_flood$Country.code[which(Hanze_flood$Year>1980)])

# matching evflood2 with the initial hanza db
Hanze_flood_1980=Hanze_flood[which(Hanze_flood$Year>1980),]
Hanze_flood_impact=full_join(Hanze_flood_1980,Ev_flood3,by=c("ï..ID"="EID"))
Hanze_flood_impact$precon="NO"
Hanze_flood_impact$precon[which(Hanze_flood_impact$NpreEv>0)]="YES"
Hanze_flood_impact$NpreEv[which(is.na(Hanze_flood_impact$NpreEv))]=0


#Amount of flash vs riverine floods detected
Hanze_flash=Hanze_flood_impact[which(Hanze_flood_impact$Type=="Flash"),]
length(which(Hanze_flash$NpreEv>0))/length(Hanze_flash$Country.code)

Hanze_riverine=Hanze_flood_impact[which(Hanze_flood_impact$Type=="River"),]
length(which(Hanze_riverine$NpreEv>0))/length(Hanze_riverine$Country.code)



#remove coastal floods
#Hanze_floodNC=Hanze_flood[-which(Hanze_flood$Type=="Coastal"),]
Hanze_flood_impactM=Hanze_flood_impact[which(!is.na(Hanze_flood_impact$Losses..2020.euro.)),]
#first boxplot
meds <- c(by(Hanze_flood_impactM$Area.flooded, Hanze_flood_impactM$precon, length))
med2 <- c(by(Hanze_flood_impact$Losses..2020.euro., Hanze_flood_impact$precon, median,na.rm=T))
#q <- c(by(ValidSY$skill, ValidSY$UpAgroup, quantile))
merdecol=match(Hanze_flood_impact$precon,names(meds))
Hanze_flood_impact$col=merdecol
lab1=c("yes","no")
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))

colo=c("1"="lightblue","2"="darkblue")

p1<-ggplot(Hanze_flood_impact, aes(x=factor(precon), y=Losses..2020.euro.)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=factor(col)),linewidth=0.8,outlier.alpha = 0.4)+
  # scale_y_continuous(limits = c(-0.5,1),name="KGE'",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  # scale_x_discrete(labels=c("1" = "100-200", "2" = "200-500",
  #                           "3" = "500-1000","4" = "1000-10 000",
  #                           "5" = "10 000-100 000","6" = ">100 000"),name="Catchment Area (km2)")+
  scale_y_log10(name = "Losses (2020 euro milions, log scale)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),  # Set the log breaks
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +  # Use log scale labels
  scale_fill_manual(values = colo, name="preconditioned flood", labels=lab1) +
  xlab("High flow anomaly in the 90 days prior to flood event")+
  
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  
  geom_text(data=data.frame(), aes(x=names(meds), y=med2+0.5*med2, label=round(meds)), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p1

#second boxplot by amound of previous high flow anomalies in previous 3 months

#3 and above are the same category
Hanze_flood_impact$preflood=round(Hanze_flood_impact$NpreEv)
Hanze_flood_impact$preflood[which(Hanze_flood_impact$preflood>3)]=3
Hanze_floodNC=Hanze_flood_impact
Hanze_flood_impactM=Hanze_flood_impact[which(!is.na(Hanze_flood_impact$Losses..2020.euro.)),]
meds <- c(by(Hanze_flood_impactM$Area.flooded, Hanze_flood_impactM$preflood, length))
med2 <- c(by(Hanze_floodNC$Losses..2020.euro., Hanze_floodNC$preflood, median,na.rm=T))
#q <- c(by(ValidSY$skill, ValidSY$UpAgroup, quantile))
merdecol=match(Hanze_floodNC$preflood,names(meds))
Hanze_floodNC$col=merdecol
lab1=c("0","1","2","3+")
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))

mean_val <- c(by(Hanze_floodNC$Losses..2020.euro., Hanze_floodNC$preflood, mean,na.rm=T))
np=names(mean_val)
mean_val=as.data.frame(mean_val)
mean_val$preflood=np
names(mean_val)[1]="mean_fat"

p1<-ggplot(Hanze_floodNC, aes(x=factor(preflood), y=Losses..2020.euro.)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.4)+
  # scale_y_continuous(limits = c(-0.5,1),name="KGE'",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  geom_point(data = mean_val, aes(x = preflood, y = mean_fat), 
             position = position_dodge(.9), size = 4, shape = 21, fill = "red", color = "white") +
  scale_x_discrete(labels=lab1,name="Number of preconditioning high flows anomlies (90 days)")+
  scale_y_log10(name = "Flood Losses (2020 euro, log scale)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),  # Set the log breaks
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +  # Use log scale labels
  #scale_fill_manual(values = colo, name="preconditioned flood", labels=lab1) +
  
  scale_fill_gradientn(
    colors=palet, n.breaks=4) +

  geom_text(data=data.frame(), aes(x=names(meds), y=med2+5e7, label=meds), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p1

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/bxplot_highflows.jpg"), width=30, height=20, units=c("cm"),dpi=600) 

#fatalities
med3 <- c(by(Hanze_floodNC$Fatalities+1, Hanze_floodNC$preflood, median,na.rm=T))
mean_val <- c(by(Hanze_floodNC$Fatalities+1, Hanze_floodNC$preflood, mean,na.rm=T))
np=names(mean_val)
mean_val=as.data.frame(mean_val)
mean_val$preflood=np
names(mean_val)[1]="mean_fat"
p2<-ggplot(Hanze_floodNC, aes(x=factor(preflood), y=Fatalities+1)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.4)+
  # Add points representing the mean for each group
  geom_point(data = mean_val, aes(x = preflood, y = mean_fat), 
             position = position_dodge(.9), size = 4, shape = 21, fill = "red", color = "white") +
  # scale_y_continuous(limits = c(-0.5,1),name="KGE'",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  scale_x_discrete(labels=lab1,name="Number of preconditioning high flow anomlies (90 days)")+
  # scale_y_continuous(name = "Fatalities", 
  #               breaks = c(10,50,100,200),  # Set the log breaks
  #               labels = c(10,50,100,200)) +  # Use log scale labels
  scale_y_log10(name = "Fatalities", 
                breaks = c(1,2,6,11,101),  # Set the log breaks
                labels = c(0,1,5,10,100)) +  # Use log scale labels
  #scale_fill_manual(values = colo, name="preconditioned flood", labels=lab1) +
  
  scale_fill_gradientn(
    colors=palet, n.breaks=4) +
  
  geom_text(data=data.frame(), aes(x=names(meds), y=med3+.5, label=meds), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p2

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/bxplot_highflows_fatal.jpg"),p2, width=30, height=20, units=c("cm"),dpi=600) 

plot(Hanze_floodNC$preflood,Hanze_floodNC$Fatalities,log="y")




#Number of events per year

#aggregate single and compound floods
Flood_yagg=aggregate(list(losses=Hanze_floodNC$Losses..2020.euro.),
                     by = list(year=Hanze_floodNC$Year,precon=Hanze_floodNC$precon),
                     FUN = function(x) c(l=length(x),s=sum(x,na.rm=T)))
Flood_yagg <- do.call(data.frame, Flood_yagg)


#new variable for the ratio
vrat=data.frame(year=Flood_yagg$year[which(Flood_yagg$precon=="YES")],
                cflood=Flood_yagg$losses.s[which(Flood_yagg$precon=="YES")],
                sflood=Flood_yagg$losses.s[which(Flood_yagg$precon=="NO")])
vrat$ratio=vrat$cflood/(vrat$cflood+vrat$sflood)*100

plot(vrat$ratio,type="l")

colorn = c("YES" ='darkblue',"NO" ='skyblue')
xlabs=seq(1950,2020,10)
clabels=c("Single","Precondition")
ggplot() +
  geom_bar(data=Flood_yagg, aes(x = year, y = losses.l, fill = factor(precon)),
           position="stack", stat="identity",alpha=0.6) +
  scale_fill_manual(values = colorn, name="Flood events",labels=clabels) +
  scale_y_continuous(
    name = "Count of Flood events",
    breaks=seq(0,100,20),
    minor_breaks = seq(-100,100,10),
    sec.axis = sec_axis( transform=~.*1, name="Percentage of losses from precondionned floods",
                         breaks=seq(0,100,10))
  )+
  geom_line(data=vrat, aes(x=year, y = 1*(ratio)), color="black",lwd=1.5) +
  # guides(color = guide_legend(override.aes = list(color = colorz)))+
  scale_x_continuous(breaks=xlabs,labels=xlabs,name="Years")+
  # geom_linerange(data=trtF,aes(x=year, ymin=fac*(cq1-1e-5),ymax=fac*cq2,color = factor(driver),group=factor(driver)),
  #                position = position_dodge2(width = 8),lwd=1) +
  #scale_color_manual(values = colorz, name="Drivers",labels=clabels) +
  theme(axis.title=element_text(size=14, face="bold"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(size=12,face="bold"),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        axis.ticks.y = element_blank(),  # Remove x-axis tickmarks
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_line(color = "lightgray"),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/values_in_time.jpg"), width=30, height=20, units=c("cm"),dpi=600) 




#Only for flash floods

Hanze_flash=Hanze_floodNC[which(Hanze_floodNC$Type=="River"),]
Hanze_flash_onlyLoss=Hanze_flash[-which(is.na(Hanze_flash$Losses..2020.euro.)),]
meds <- c(by(Hanze_flash_onlyLoss$Area.flooded, Hanze_flash_onlyLoss$preflood, length))
med2 <- c(by(Hanze_flash$Losses..2020.euro., Hanze_flash$preflood, median,na.rm=T))
#q <- c(by(ValidSY$skill, ValidSY$UpAgroup, quantile))
merdecol=match(Hanze_flash$preflood,names(meds))
Hanze_flash$col=merdecol
lab1=c("0","1","2","3+")
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))

mean_val <- c(by(Hanze_flash$Losses..2020.euro., Hanze_flash$preflood, mean,na.rm=T))
np=names(mean_val)
mean_val=as.data.frame(mean_val)
mean_val$preflood=np
names(mean_val)[1]="mean_fat"

p1<-ggplot(Hanze_flash, aes(x=factor(preflood), y=Losses..2020.euro.)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.4)+
  # scale_y_continuous(limits = c(-0.5,1),name="KGE'",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  geom_point(data = mean_val, aes(x = preflood, y = mean_fat), 
             position = position_dodge(.9), size = 4, shape = 21, fill = "red", color = "white") +
  scale_x_discrete(labels=lab1,name="Number of preconditioning high flows anomalies (90 days)")+
  scale_y_log10(name = "Flash Flood Losses (2020 euro, log scale)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),  # Set the log breaks
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +  # Use log scale labels
  #scale_fill_manual(values = colo, name="preconditioned flood", labels=lab1) +
  
  scale_fill_gradientn(
    colors=palet, n.breaks=4) +
  
  geom_text(data=data.frame(), aes(x=names(meds), y=med2+5e7, label=meds), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p1


#######  same for low flow precondition #########


#aggregate by event ID
Ev_lowfloods=aggregate(list(time=Events_precoflood$time),
                    by = list(NID=Events_precoflood$NUTID,EID=Events_precoflood$eventID),
                    FUN = function(x) c(l=length(x)))
names(Ev_lowfloods)[3]="N_PreEv"

Ev_lowfloods$NEID=paste0(Ev_lowfloods$EID,Ev_lowfloods$NID)

Hanze_flood_BN$NEID=paste0(Hanze_flood_BN$ï..ID,Hanze_flood_BN$NUTl)

Ev_lowflood2=full_join(Ev_floods,Hanze_flood_BN,by="NEID")

Ev_lowflood2$N_PreEv[which(is.na(Ev_lowflood2$N_PreEv))]=0
Ev_lowflood2$NID[which(is.na(Ev_lowflood2$NID))]=Ev_lowflood2$NUTl[which(is.na(Ev_lowflood2$EID))]
Ev_lowflood2$EID[which(is.na(Ev_lowflood2$EID))]=Ev_lowflood2$ï..ID[which(is.na(Ev_lowflood2$EID))]

#Identify events for which the previous event is the event itself

Ev_lowflood2$PC=Ev_lowflood2$N_PreEv
Ev_lowflood2$PC[which(Ev_lowflood2$PC>0)]=1

Regio_mplot=aggregate(list(time=Ev_lowflood2$Losses..2020.euro.),
                      by = list(NID=Ev_lowflood2$NID,Preflood=Ev_lowflood2$PC),
                      FUN = function(x) c(l=length(x),mean=mean(x,na.rm=T),median=median(x,na.rm=T)))
Regio_mplot <- do.call(data.frame, Regio_mplot)

Regio_mplot0=Regio_mplot[which(Regio_mplot$Preflood==0),]
Regio_mplot1=Regio_mplot[which(Regio_mplot$Preflood==1),]
Regio_mplot0=full_join(Regio_mplot1,Regio_mplot0,by="NID")
Regio_mplot0$Ratio=Regio_mplot0$time.mean.x/Regio_mplot0$time.mean.y
#NUTS3 maps
Regio_mplot=inner_join(Regio_mplot0,NUTS3,by=c("NID"="NUTS_ID"))



fdp=ggplot(basemap) +
  geom_sf(fill="gray95")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=Regio_mplot,aes(fill=time.l.x,geometry=geometry),color="gray12",alpha=1,size=0.25,stroke=0,shape=15)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet,
    oob = scales::squish,na.value="gray30", name="Loss ratio") +
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("Ratio of losses from preconditionned flood events over single flood events (low flows)")
fdp

#Nplot <- st_transform(Regio_mplot, crs = 3035)


nuts1 <- st_read("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/NUTS_RG_20M_2021_3035.shp")


Regio_mplot=st_as_sf(Regio_mplot)
nuts1=nuts1[which(nuts1$LEVL_CODE==1),]

n3to1 <- st_join( nuts1, Regio_mplot,join = st_intersects)


#aggregate byy nuts1


Regio_n1plot_cp=aggregate(list(val=n3to1$time.mean.x),
                       by = list(NID=n3to1$NUTS_ID),
                       FUN = function(x) c(l=length(x),sum=sum(x,na.rm=T),mean=mean(x,na.rm=T),median=median(x,na.rm=T)))
Regio_n1plot_cp <- do.call(data.frame, Regio_n1plot_cp)


Regio_n1plot_sg=aggregate(list(val=n3to1$time.mean.y),
                          by = list(NID=n3to1$NUTS_ID),
                          FUN = function(x) c(l=length(x),sum=sum(x,na.rm=T),mean=mean(x,na.rm=T),median=median(x,na.rm=T)))
Regio_n1plot_sg <- do.call(data.frame, Regio_n1plot_sg)

Regio_n1plot_cp$ratio=(Regio_n1plot_cp$val.sum+1)/(Regio_n1plot_sg$val.sum+1)
Regio_n1plot_cp$ratio[which(Regio_n1plot_cp$val.sum==0 & Regio_n1plot_sg$val.sum==0)]=NA
Regio_n1plot_cp$sum = (Regio_n1plot_cp$val.sum+1)+(Regio_n1plot_sg$val.sum+1)
Regio_n1plot_cp$sum[which(Regio_n1plot_cp$sum==2)]=NA

Regio_n1plot=inner_join(Regio_n1plot_cp,nuts1,by=c("NID"="NUTS_ID"))






palet=c(hcl.colors(11, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
paletf=c(hcl.colors(11, palette = "RdBu", alpha = NULL, rev = T, fixup = TRUE))

pl1=ggplot(basemap) +
  geom_sf(fill="gray95")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=Regio_n1plot,aes(fill=ratio,geometry=geometry),color="gray12",alpha=1,size=0.25,stroke=0,shape=15)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=paletf,limits=c(1e-1,1e1),trans="log", breaks=c(1e-2,1e-1,1,10,100),
    oob = scales::squish,na.value="gray30", name="Loss ratio") +
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("Ratio of losses from preconditionned flood events over single flood events (low flows)")
pl1


ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/LossRatioMap_lowflow.jpg"),pl1, width=30, height=25, units=c("cm"),dpi=800) 

pl2=ggplot(basemap) +
  geom_sf(fill="gray95")+
  geom_sf(fill=NA, color="grey") +
  geom_sf(data=Regio_n1plot,aes(fill=sum/1000000,geometry=geometry),color="gray12",alpha=1,size=0.25,stroke=0,shape=15)+ 
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet,trans="log", breaks=c(1,10,100,1000,10000,100000),
    oob = scales::squish,na.value="gray30", name="Total losses (million Euros)") +
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("Total losses")
pl2









Hanze_floodEV=aggregate(list(time=Hanze_flood_BN$Year),
                        by = list(NID=Hanze_flood_BN$NUTl),
                        FUN = function(x) c(l=length(x)))



Ev_lowflood3=aggregate(list(Ev_lowflood2$N_PreEv),
                    by = list(EID=Ev_lowflood2$EID),
                    FUN = function(x) c(l=length(x),s=mean(x)))
Ev_lowflood3 <- do.call(data.frame, Ev_lowflood3)
names(Ev_lowflood3)[c(2,3)]=c("Nregion","NpreEv")

length(Ev_lowflood3$EID[which(Ev_lowflood3$NpreEv==0)])/length(Hanze_flood$Country.code[which(Hanze_flood$Year>1980)])

# matching evflood2 with the initial hanza db
Hanze_flood_1980=Hanze_flood[which(Hanze_flood$Year>1980),]
Hanze_flood_impact=full_join(Hanze_flood_1980,Ev_lowflood3,by=c("ï..ID"="EID"))
Hanze_flood_impact$precon="NO"
Hanze_flood_impact$precon[which(Hanze_flood_impact$NpreEv>0)]="YES"
Hanze_flood_impact$NpreEv[which(is.na(Hanze_flood_impact$NpreEv))]=0


#Amount of flash vs riverine floods detected
Hanze_flash=Hanze_flood_impact[which(Hanze_flood_impact$Type=="Flash"),]
length(which(Hanze_flash$NpreEv>0))/length(Hanze_flash$Country.code)

Hanze_riverine=Hanze_flood_impact[which(Hanze_flood_impact$Type=="River"),]
length(which(Hanze_riverine$NpreEv>0))/length(Hanze_riverine$Country.code)



#remove coastal floods
#Hanze_floodNC=Hanze_flood[-which(Hanze_flood$Type=="Coastal"),]
Hanze_flood_impactM=Hanze_flood_impact[which(!is.na(Hanze_flood_impact$Losses..2020.euro.)),]
#first boxplot
meds <- c(by(Hanze_flood_impactM$Area.flooded, Hanze_flood_impactM$precon, length))
med2 <- c(by(Hanze_flood_impact$Losses..2020.euro., Hanze_flood_impact$precon, median,na.rm=T))
#q <- c(by(ValidSY$skill, ValidSY$UpAgroup, quantile))
merdecol=match(Hanze_flood_impact$precon,names(meds))
Hanze_flood_impact$col=merdecol
lab1=c("yes","no")
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))

colo=c("1"="lightblue","2"="darkblue")

p1<-ggplot(Hanze_flood_impact, aes(x=factor(precon), y=Losses..2020.euro.)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=factor(col)),linewidth=0.8,outlier.alpha = 0.4)+
  # scale_y_continuous(limits = c(-0.5,1),name="KGE'",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  # scale_x_discrete(labels=c("1" = "100-200", "2" = "200-500",
  #                           "3" = "500-1000","4" = "1000-10 000",
  #                           "5" = "10 000-100 000","6" = ">100 000"),name="Catchment Area (km2)")+
  scale_y_log10(name = "Losses (2020 euro milions, log scale)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),  # Set the log breaks
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +  # Use log scale labels
  scale_fill_manual(values = colo, name="preconditioned flood", labels=lab1) +
  xlab("High flow anomaly in the 90 days prior to flood event")+
  
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  
  geom_text(data=data.frame(), aes(x=names(meds), y=med2+0.5*med2, label=round(meds)), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p1

#second boxplot by amound of previous high flow anomalies in previous 3 months

#3 and above are the same category
Hanze_flood_impact$preflood=round(Hanze_flood_impact$NpreEv)
Hanze_flood_impact$preflood[which(Hanze_flood_impact$preflood>3)]=3
Hanze_floodNC=Hanze_flood_impact
Hanze_flood_impactM=Hanze_flood_impact[which(!is.na(Hanze_flood_impact$Losses..2020.euro.)),]
meds <- c(by(Hanze_flood_impactM$Area.flooded, Hanze_flood_impactM$preflood, length))
med2 <- c(by(Hanze_floodNC$Losses..2020.euro., Hanze_floodNC$preflood, median,na.rm=T))
#q <- c(by(ValidSY$skill, ValidSY$UpAgroup, quantile))
merdecol=match(Hanze_floodNC$preflood,names(meds))
Hanze_floodNC$col=merdecol
lab1=c("0","1","2","3+")
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))

mean_val <- c(by(Hanze_floodNC$Losses..2020.euro., Hanze_floodNC$preflood, mean,na.rm=T))
np=names(mean_val)
mean_val=as.data.frame(mean_val)
mean_val$preflood=np
names(mean_val)[1]="mean_fat"

p1<-ggplot(Hanze_floodNC, aes(x=factor(preflood), y=Losses..2020.euro.)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.4)+
  # scale_y_continuous(limits = c(-0.5,1),name="KGE'",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  geom_point(data = mean_val, aes(x = preflood, y = mean_fat), 
             position = position_dodge(.9), size = 4, shape = 21, fill = "red", color = "white") +
  scale_x_discrete(labels=lab1,name="Number of preconditioning low flows anomlies (90 days)")+
  scale_y_log10(name = "Flood Losses (2020 euro, log scale)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),  # Set the log breaks
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +  # Use log scale labels
  #scale_fill_manual(values = colo, name="preconditioned flood", labels=lab1) +
  
  scale_fill_gradientn(
    colors=palet, n.breaks=4) +
  
  geom_text(data=data.frame(), aes(x=names(meds), y=med2+5e7, label=meds), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p1

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/bxplot_lowflows.jpg"), width=30, height=20, units=c("cm"),dpi=600) 


#fatalities
med3 <- c(by(Hanze_floodNC$Fatalities+1, Hanze_floodNC$predrought, median,na.rm=T))
mean_val <- c(by(Hanze_floodNC$Fatalities+1, Hanze_floodNC$predrought, mean,na.rm=T))
np=names(mean_val)
mean_val=as.data.frame(mean_val)
mean_val$predrought=np
names(mean_val)[1]="mean_fat"
p2<-ggplot(Hanze_floodNC, aes(x=factor(predrought), y=Fatalities+1)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.4)+
  # Add points representing the mean for each group
  geom_point(data = mean_val, aes(x = predrought, y = mean_fat), 
             position = position_dodge(.9), size = 4, shape = 21, fill = "red", color = "white") +
  # scale_y_continuous(limits = c(-0.5,1),name="KGE'",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  scale_x_discrete(labels=lab1,name="Number of preconditioning low flow anomalies (90 days)")+
  # scale_y_continuous(name = "Fatalities", 
  #               breaks = c(10,50,100,200),  # Set the log breaks
  #               labels = c(10,50,100,200)) +  # Use log scale labels
  scale_y_log10(name = "Fatalities", 
                breaks = c(1,2,6,11,101),  # Set the log breaks
                labels = c(0,1,5,10,100)) +  # Use log scale labels
  #scale_fill_manual(values = colo, name="preconditioned flood", labels=lab1) +
  
  scale_fill_gradientn(
    colors=palet, n.breaks=4) +
  
  geom_text(data=data.frame(), aes(x=names(meds), y=med3+.5, label=meds), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p2



#Only for flash floods

Hanze_flash=Hanze_floodNC[which(Hanze_floodNC$Type=="Flash"),]
Hanze_flash$predrought[which(Hanze_flash$predrought>3)]=3
Hanze_flash_onlyLoss=Hanze_flash[-which(is.na(Hanze_flash$Losses..2020.euro.)),]
meds <- c(by(Hanze_flash_onlyLoss$Area.flooded, Hanze_flash_onlyLoss$predrought, length))
med2 <- c(by(Hanze_flash$Losses..2020.euro., Hanze_flash$predrought, median,na.rm=T))
#q <- c(by(ValidSY$skill, ValidSY$UpAgroup, quantile))
merdecol=match(Hanze_flash$predrought,names(meds))
Hanze_flash$col=merdecol
lab1=c("0","1","2","3+")
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))

mean_val <- c(by(Hanze_flash$Losses..2020.euro., Hanze_flash$predrought, mean,na.rm=T))
np=names(mean_val)
mean_val=as.data.frame(mean_val)
mean_val$preflood=np
names(mean_val)[1]="mean_fat"

p1<-ggplot(Hanze_flash, aes(x=factor(predrought), y=Losses..2020.euro.)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.4)+
  # scale_y_continuous(limits = c(-0.5,1),name="KGE'",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  geom_point(data = mean_val, aes(x = preflood, y = mean_fat), 
             position = position_dodge(.9), size = 4, shape = 21, fill = "red", color = "white") +
  scale_x_discrete(labels=lab1,name="Number of preconditioning low flows anomalies (90 days)")+
  scale_y_log10(name = "Flash Flood Losses (2020 euro, log scale)", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),  # Set the log breaks
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +  # Use log scale labels
  #scale_fill_manual(values = colo, name="preconditioned flood", labels=lab1) +
  
  scale_fill_gradientn(
    colors=palet, n.breaks=4) +
  
  geom_text(data=data.frame(), aes(x=names(meds), y=med2+5e7, label=meds), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p1





#Number of events per year

#aggregate single and compound floods
Flood_yagg=aggregate(list(losses=Hanze_floodNC$Losses..2020.euro.),
                    by = list(year=Hanze_floodNC$Year,precon=Hanze_floodNC$precon),
                    FUN = function(x) c(l=length(x),s=sum(x,na.rm=T)))
Flood_yagg <- do.call(data.frame, Flood_yagg)


#new variable for the ratio
vrat=data.frame(year=Flood_yagg$year[which(Flood_yagg$precon=="YES")],
                cflood=Flood_yagg$losses.s[which(Flood_yagg$precon=="YES")],
                sflood=Flood_yagg$losses.s[which(Flood_yagg$precon=="NO")])
vrat$ratio=vrat$cflood/(vrat$cflood+vrat$sflood)*100

plot(vrat$ratio,type="l")

colorn = c("YES" ='darkblue',"NO" ='skyblue')
xlabs=seq(1950,2020,10)
clabels=c("Single","Precondition")
ggplot() +
  geom_bar(data=Flood_yagg, aes(x = year, y = losses.l, fill = factor(precon)),
           position="stack", stat="identity",alpha=0.6) +
  scale_fill_manual(values = colorn, name="Flood events",labels=clabels) +
  scale_y_continuous(
    name = "Count of Flood events",
    breaks=seq(0,100,20),
    minor_breaks = seq(-100,100,10),
    sec.axis = sec_axis( transform=~.*2, name="Percentage of losses from precondionned floods",
                         breaks=seq(0,100,10))
  )+
  geom_line(data=vrat, aes(x=year, y = 0.5*(ratio)), color="black",lwd=1.5) +
  # guides(color = guide_legend(override.aes = list(color = colorz)))+
  scale_x_continuous(breaks=xlabs,labels=xlabs,name="Years")+
  # geom_linerange(data=trtF,aes(x=year, ymin=fac*(cq1-1e-5),ymax=fac*cq2,color = factor(driver),group=factor(driver)),
  #                position = position_dodge2(width = 8),lwd=1) +
  #scale_color_manual(values = colorz, name="Drivers",labels=clabels) +
  theme(axis.title=element_text(size=14, face="bold"),
        axis.text = element_text(size=12),
        axis.text.x = element_text(size=12,face="bold"),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        axis.ticks.y = element_blank(),  # Remove x-axis tickmarks
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_line(color = "lightgray"),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))






#Aggregating number of compound and not compound events per year







Events_prefflood=Events_prefflood[unique(Events_prefflood$eventID),]

#I have to create a new table with NID, EiD and losses



boxplot(Events_prefflood$Losses..2020.euro.,log="y")

boxplot(Hanze_flood$Losses..2020.euro.,log="y")

#I need to go back to the inital HANZE setting


df_POTh_large=c()
df_POTl_large=c()
Hanze_flood_BN=c()
High_events_df=c()
Low_events_df=c()
for (id in 1:length(Hanze_ids)){
  print(id)
  NUTl=NutVector[id]
  myloc=NR_Qd_total[,id]
  if (length(which(!is.na(myloc)))>0){
    ms=data.frame(time=NRtime,data=myloc)
    
    POTdata_high <- get_POTdata_high(95,ms)
    POTdata_high$enddate=as.POSIXct(POTdata_high$endpeaks*24*3600-3600,origin="1980-01-01 00:00:00")
    POTdata_low <- get_POTdata_sliding_low(05, ms)
    
    #saving raw events detected
    Low_events=POTdata_low$low_peaks_df
    Low_events$NUT=NUTl
    High_events=POTdata_high
    High_events$NUT=NUTl
    
    Low_events_df=rbind(Low_events_df,Low_events)
    High_events_df=rbind(High_events_df,High_events)
    
    
    matches <- lapply(list_of_word_vectors, function(vector) NUTl %in% vector)
    tm=which(matches==T)
    Hanze_Nut=Hanze_flood[tm,]
    Hanze_Nut=Hanze_Nut[which(Hanze_Nut$Year>1980),]
    if(length(Hanze_Nut$Country.code)>0){
      Hanze_bNut=data.frame(NUTl,Hanze_Nut)
      Hanze_flood_BN=rbind(Hanze_flood_BN,Hanze_bNut)
    }
    #now check overlap with mt detected anomalies
    df_POTh_keep=c()
    df_POTl_keep=c()
    for (e in 1:length(Hanze_Nut$ï..ID)){
      d1=(as.Date(POTdata_high$enddate)-as.Date(Hanze_Nut$End.date[e]))
      keep=which(d1<7 & d1>-90)
      if (length(keep)>0){
        POTh_keep=POTdata_high[keep,]
        POTh_keep$dtime=d1[keep]
        POTh_keep$NUTID=NUTl
        POTh_keep$eventID=Hanze_Nut$ï..ID[e]
        df_POTh_keep=rbind(df_POTh_keep,POTh_keep)
      } else{
        
      }
      
      
      d2=(as.Date(POTdata_low$low_peaks_df$endpeaks)-as.Date(Hanze_Nut$End.date[e]))
      keep=which(d2<0 & d2>-90)
      
      if (length(keep)>0){
        POTl_keep=POTdata_low$low_peaks_df[keep,]
        POTl_keep$dtime=d2[keep]
        POTl_keep$NUTID=NUTl
        POTl_keep$eventID=Hanze_Nut$ï..ID[e]
        df_POTl_keep=rbind(df_POTl_keep,POTl_keep)
      }
    }
    
    df_POTh_keep=as.data.frame(df_POTh_keep)
    df_POTl_keep=as.data.frame(df_POTl_keep)
    
    df_POTh_large=rbind(df_POTh_large,df_POTh_keep)
    df_POTl_large=rbind(df_POTl_large,df_POTl_keep)
  }
  
}

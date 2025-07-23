#Library importation
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
library(ncdf4)
library(foreach)
library(doParallel)
library(raster)
library(sf)
library(exactextractr)

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

#Functions
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



#Functions declaration -----

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
  
  # Create the data frame directly with the extracted data
  POTdata_df <- data.frame(
    threshold = thrsd,         # Same threshold value for all rows
    percentile = pcts,         # Same percentile value for all rows
    peaks = pks,               # Peak values
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
      
      POTdata_df <- rbind(POTdata_df, data.frame(
        threshold = thresholds[doy[ipeaks]],
        percentile = pcts,
        low_peaks = ms[ipeaks, 2],  # The lowest point in the cluster
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

#Season computation functions
seasony=function(x){
  theta=x*(2*pi/365.25)
  #plot(cos(theta),sin(theta),xlim=c(-1,1),ylim=c(-1,1))
  xi=1/(length(theta))*sum(cos(na.omit(theta)))
  yi=1/(length(theta))*sum(sin(na.omit(theta)))
  if (xi<=0){
    Di=(atan(yi/xi)+pi)*(365.25/(2*pi))
  }else if(xi>0 & yi>=0){
    Di=(atan(yi/xi))*(365.25/(2*pi))
  }else if(xi>0 & yi<0){
    Di=(atan(yi/xi)+2*pi)*(365.25/(2*pi))
  }
  R=sqrt(xi^2+yi^2)
  return(c(Di,R))
}
season1=function(x){
  l1=length(which(!is.na(x)))
  if(l1>0){
    x=x[which(!is.na(x))]
    theta=x*(2*pi/365.25)
    # plot(theta)
    
    xi=1/(length(theta))*sum(cos(theta))
    yi=1/(length(theta))*sum(sin(theta))
    if (xi<=0){
      Di=(atan(yi/xi)+pi)*(365.25/(2*pi))
    }else if(xi>0 & yi>=0){
      Di=(atan(yi/xi))*(365.25/(2*pi))
    }else if(xi>0 & yi<0){
      Di=(atan(yi/xi)+2*pi)*(365.25/(2*pi))
    }
    R=sqrt(xi^2+yi^2)
  }else{Di=NA}
  return(Di)
}
season2=function(x){
  l1=length(which(!is.na(x)))
  if(l1>0){
    x=x[which(!is.na(x))]
    theta=x*(2*pi/365.25)
    # plot(theta)
    
    xi=1/(length(theta))*sum(cos(theta))
    yi=1/(length(theta))*sum(sin(theta))
    if (xi<=0){
      Di=(atan(yi/xi)+pi)*(365.25/(2*pi))
    }else if(xi>0 & yi>=0){
      Di=(atan(yi/xi))*(365.25/(2*pi))
    }else if(xi>0 & yi<0){
      Di=(atan(yi/xi)+2*pi)*(365.25/(2*pi))
    }
    R=sqrt(xi^2+yi^2)
  }else{Di=NA
  R=NA}
  return(R)
}

direction_labeller <- function(x){
  ifelse(x %% 30.5 == 0, c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct',"Nov",'Dec')[1+(as.integer(x/30.5) %% 12)], '')
}

elbow_finder <- function(x_values, y_values) {
  # Max values to create line
  max_x_x <- min(x_values)
  max_x_y <- y_values[which.min(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))
  
  # Creating straight line between the max values
  sl=approx(max_df$y, max_df$x,n=15)
  
  xcurve=sl$y/max(sl$y)
  ycurve=sl$x/max(sl$x)
  x_values2=x_values/max(x_values)
  y_values2=y_values/max(y_values)
  plot(xcurve,ycurve)
  lines(x_values2,y_values2,col=2)
  # Distance from point to line
  x_min_dist<-c()
  y_min_dist<-c()
  distances <- c()
  for(i in 1:length(x_values)) {
    distancex<-c()
    for(j in 1:length(x_values)){
      #distancex <- c(distancex, abs(coef(fit)[2]*x_values[j] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
      print(sqrt((x_values2[i]-xcurve[j])^2+(y_values2[i]-ycurve[j])^2))
      distancex<-c(distancex,sqrt((x_values2[i]-xcurve[j])^2+(y_values2[i]-ycurve[j])^2))
    }
    plot(distancex)
    distances<-c(distances,min(distancex, na.rm=T))
    print(which.min(distancex))
    x_min_dist <- c(x_min_dist,xcurve[which.min(distancex)])
    y_min_dist <- c(y_min_dist,ycurve[which.min(distancex)])
  }
  
  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]
  xcmax <- x_min_dist[which.max(distances)]
  ycmax <- y_min_dist[which.max(distances)]
  
  
  poncu<-c()
  plot(x_values,distances)
  
  return(c(x_max_dist, y_max_dist,xcmax,ycmax))
}

# Function to calculate proportions above each percentile
calculate_proportions <- function(df, percentiles, damage_percentiles) {
  proportions <- c()
  
  total_loss=sum(df$damage)
  for (i in seq_along(percentiles)) {
    threshold <- damage_percentiles[i]
    total_events <- sum(df$damage >= threshold)
    sum_events <- sum(df$damage[which(df$damage >= threshold)])/total_loss
    #lenght_percentile<- c(length_percentile,length_events)
    for (cat in unique(df$category)) {
      above_threshold <- sum(df$category == cat & df$damage >= threshold)
      prop <- above_threshold / total_events
      proportions <- rbind(proportions, data.frame(percentile = percentiles[i], category = cat, proportion = prop,
                                                   length=total_events, damage_perc=sum_events))
    }
  }
  
  return(proportions)
}


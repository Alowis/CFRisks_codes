

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
  
  }
}

#save outputs from here
save(High_events_df,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/events_highflows_1980-2020_v3.Rdata"))
save(Low_events_df,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/events_lowflows_1980-2020_v2.Rdata"))

#new loop for floods with varying threshold
load(file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/events_highflows_1980-2020.Rdata"))
load(file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/events_lowflows_1980-2020.Rdata"))


#save events in csv
write.csv(High_events_df,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Hflow_Q_1980-2020.csv"))
write.csv(Low_events_df,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Lflow_Q_1980-2020.csv"))

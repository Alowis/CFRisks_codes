
source("D:/tilloal/Documents/01_Projects/RiskDynamics/CFRisk_codes/functions_CFRisks.R")
#Set data directory
hydroDir<-("D:/tilloal/Documents/01_Projects/RiskDynamics/CFRisks_data/")

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
    outhloc=outhybas
    outf=rbind(outf,outhloc)
  }
}


### NUTS3 ----
NUTS3 <- read_sf(dsn = paste0(hydroDir,"NUTS3/NUTS3_modified.shp"))
GridNUTS3=raster( paste0(hydroDir,"NUTS3/NUTS3_Raster3ID.tif"))
GN3=as.data.frame(GridNUTS3,xy=T)
GN3=GN3[which(!is.na(GN3[,3])),]
GN3$llcoord=paste(round(GN3$x,4),round(GN3$y,4),sep=" ")
GN3_riv=right_join(GN3,outf,by= c("llcoord"="latlong"))

GridNUTS2=raster( paste0(hydroDir,"NUTS3/NUTS3_Raster2ID.tif"))
GN2=as.data.frame(GridNUTS2,xy=T)
GN2=GN2[which(!is.na(GN2[,3])),]
GN2$llcoord=paste(round(GN2$x,4),round(GN2$y,4),sep=" ")
GN2_riv=right_join(GN2,outf,by= c("llcoord"="latlong"))

GNF=right_join(GN3,GN2_riv,by="llcoord")
GNUTS3sf=fortify(NUTS3)
GNFx=GNF[which(is.na(GNF$NUTS3_Raster3ID)),]


#load UpArea
#load upstream area
outletname="upArea_European_01min.nc"
dir=hydroDir
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(hydroDir,outletname,outf)
head(UpArea)


#Parallel computing for aggregation ----

# Register parallel backend to use multiple cores
no_cores <- detectCores() - 1 # Leave one core free
cl <- makeCluster(no_cores)
registerDoParallel(cl)

yseq=seq(1980,2020)
NUTQspT=c()
lnut=length(NUTS3$NUTS_ID)
datafolder="D:/tilloal/Documents/06_Floodrivers"

results <- foreach(year = yseq, .packages = c("dplyr","ncdf4", "raster", "exactextractr", "sf")) %dopar%
{
  
  # Open the NetCDF file
  print(year)
  nc <- nc_open(paste0(datafolder,"/HERA2/dis.HERA2_",year,".nc"))
  
  # Extract the dimensions of the file
  lon_dim <- nc$dim[["lon"]]
  lat_dim <- nc$dim[["lat"]]
  time_dim <- nc$dim[["time"]]
  name.vb=names(nc[['var']])
  namev=name.vb[2]
  time <- ncvar_get(nc,"time")
  lt=length(time)
  
  # Extract the longitude and latitude values
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  outll=expand.grid(lon,lat)
  
  first_day <- ncvar_get(nc, namev, 
                         start = c(1, 1, 1), 
                         count = c(lon_dim$len, lat_dim$len, 1))
  
  # Create a data frame from the first day
  dx <- data.frame(lon = as.vector(outll$Var1), 
                   lat = as.vector(outll$Var2), 
                   value = as.vector(first_day))
  
  dx2=dx[-which(is.na(dx$value)),]
  
  #divide by upArea to have specific discharge
  dx2$ll=paste(round(dx2$lon,4),round(dx2$lat,4))
  mup=match(UpArea$latlong,dx2$ll)
  NUTQsp=matrix(ncol=3, nrow=lt*length(NUTS3$NUTS_ID))
  for (t in 1:lt){
    print(t)
    # Extract the first day of the file
    dataday <- ncvar_get(nc, namev, 
                         start = c(1, 1, t), 
                         count = c(lon_dim$len, lat_dim$len, 1))
    
    
    # Create a data frame from the first day
    df <- data.frame(lon = as.vector(outll$Var1), 
                     lat = as.vector(outll$Var2), 
                     value = as.vector(dataday))
    
    df2=df[-which(is.na(df$value)),]
    
    #divide by upArea to have specific discharge
    df2$upa=UpArea$upa[mup]
    df2$qsp=df2$value/df2$upa*1000
    
    # Create a raster layer from the data frame
    r <- rasterFromXYZ(df2[, c("lon", "lat", "qsp")], 
                       crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    
    mean_NUTS3<- exact_extract(r, NUTS3, 'mean')
    
    #extraction of most important features:
    
    
    NUTsave=cbind(NUTS3[,9],mean_NUTS3)
    st_geometry(NUTsave)<-NULL
    NUTsave=as.data.frame(NUTsave)
    NUTsave$time=time[t]
    NUTsave=as.matrix(NUTsave)
    
    NUTQsp[c((1+(t-1)*lnut):(t*lnut)),]=NUTsave
  }
  
  file_name <- paste0("/NUTSQ/NUTS3xx_", year, ".csv")
  write.csv(NUTQsp,file=paste0(hydroDir,file_name))
  
  #NUTQspT=rbind(NUTQspT,NUTQsp)
  
}


stopCluster(cl)

#

for (year in yseq)
{
  # Open the NetCDF file
  print(year)
  nc <- nc_open(paste0(datafolder,"/HERA2/dis.HERA2_",year,".nc"))
  
  # Extract the dimensions of the file
  lon_dim <- nc$dim[["lon"]]
  lat_dim <- nc$dim[["lat"]]
  time_dim <- nc$dim[["time"]]
  name.vb=names(nc[['var']])
  namev=name.vb[2]
  time <- ncvar_get(nc,"time")
  lt=length(time)
  
  
  # Extract the longitude and latitude values
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  outll=expand.grid(lon,lat)
  
  first_day <- ncvar_get(nc, namev, 
                       start = c(1, 1, 1), 
                       count = c(lon_dim$len, lat_dim$len, 1))
  
  # Create a data frame from the first day
  dx <- data.frame(lon = as.vector(outll$Var1), 
                   lat = as.vector(outll$Var2), 
                   value = as.vector(first_day))
  
  dx2=dx[-which(is.na(dx$value)),]
  
  #divide by upArea to have specific discharge
  dx2$ll=paste(round(dx2$lon,4),round(dx2$lat,4))
  mup=match(UpArea$latlong,dx2$ll)
  NUTQsp=matrix(ncol=3, nrow=lt*length(NUTS3$NUTS_ID))
  for (t in 1:lt){
    print(t)
    # Extract the first day of the file
    dataday <- ncvar_get(nc, namev, 
                           start = c(1, 1, t), 
                           count = c(lon_dim$len, lat_dim$len, 1))
    
    
    # Create a data frame from the first day
    df <- data.frame(lon = as.vector(outll$Var1), 
                     lat = as.vector(outll$Var2), 
                     value = as.vector(dataday))
    
    df2=df[-which(is.na(df$value)),]
    
    #divide by upArea to have specific discharge
    df2$upa=UpArea$upa[mup]
    df2$qsp=df2$value/df2$upa*1000
    
    # Create a raster layer from the data frame
    r <- rasterFromXYZ(df2[, c("lon", "lat", "qsp")], 
                       crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    
    mean_NUTS3<- exact_extract(r, NUTS3, 'mean')
    
    # ggplot(basemap) +
    #   geom_sf(data=NUTS3,aes(fill=mean_NUTS3,geometry=geometry),color="transparent")+
    #   coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    #   labs(x="Longitude", y = "Latitude")+
    #   theme(axis.title=element_text(size=tsize),
    #         title = element_text(size=osize),
    #         axis.text=element_text(size=osize),
    #         panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
    #         panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
    #         legend.title = element_text(size=tsize),
    #         legend.text = element_text(size=osize),
    #         legend.position = "right",
    #         panel.grid.major = element_line(colour = "grey70"),
    #         panel.grid.minor = element_line(colour = "grey90"),
    #         legend.key = element_rect(fill = "transparent", colour = "transparent"),
    #         legend.key.size = unit(1, "cm"))
    
    
    #extraction of most important features:
    
    
    NUTsave=cbind(NUTS3[,9],mean_NUTS3)
    st_geometry(NUTsave)<-NULL
    NUTsave=as.data.frame(NUTsave)
    NUTsave$time=time[t]
    NUTsave=as.matrix(NUTsave)
    
    NUTQsp[c((1+(t-1)*lnut):(t*lnut)),]=NUTsave
  }
  file_name <- paste0("/NUTSQ/NUTS3x_", year, ".csv")
  write.csv(NUTQsp,file=paste0(hydroDir,file_name))

}



NID=NUTS3$NUTS_ID
qplot=NUTQsp[which(NUTQsp$NUTS_ID==NID[400]),]

plot(time,qplot$mean_NUTS3)
#transform df2 into a rasterlayer
# Close the NetCDF file
nc_close(nc)





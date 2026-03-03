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

#Load my shapefile on which to aggregate



### HydroRegions ----

# GridHR=raster( paste0(hydroDir,"/HydroRegions_raster_WGS84.tif"))
# GHR=as.data.frame(GridHR,xy=T)
# GHR=GHR[which(!is.na(GHR[,3])),]
# GHR$llcoord=paste(round(GHR$x,4),round(GHR$y,4),sep=" ") 
# GHR_riv=inner_join(GHR,outf,by= c("llcoord"="latlong"))
# GHshpp <- read_sf(dsn ="Z:/ClimateRun4/nahaUsers/tilloal/HydroRegions/her_all_adjusted.shp")
# HydroRsf=fortify(GHshpp) 


### NUTS3 ----


NUTS3 <- read_sf(dsn = paste0(hydroDir,"/Countries/NUTS3/NUTS3_Extended_domain.shp"))
NUTS3$N3ID=c(1:length(NUTS3$NUTS_ID))
N2ID=unique(NUTS3$NUTS2_ID)
N2IDn=c(1:length(N2ID))
mati=match(NUTS3$NUTS2_ID,N2ID)
NUTS3$N2ID=N2IDn[mati]
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



#load UpArea
#load upstream area
main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
outletname="/GIS/upArea_European_01min.nc"
dir=valid_path
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(valid_path,outletname,outf)
head(UpArea)



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


#load data to be aggregated

#Discharge data from HERA

#load yearly file

# Install and load the ncdf4 package
library(ncdf4)
library(foreach)
library(doParallel)
library(raster)
library(sf)
library(exactextractr)

# Register parallel backend to use multiple cores
no_cores <- detectCores() - 1 # Leave one core free
cl <- makeCluster(no_cores)
registerDoParallel(cl)

yseq=seq(1980,2020)
# yseq=seq(1980,1985)
NUTQspT=c()


lnut=length(NUTS3$NUTS_ID)
# 
results <- foreach(year = yseq, .packages = c("dplyr","ncdf4", "raster", "exactextractr", "sf")) %dopar%
{
  
  # Open the NetCDF file
  print(year)
  nc <- nc_open(paste0("D:/tilloal/Documents/06_Floodrivers/HERA/dis.HERA_",year,".nc"))
  
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
  
  file_name <- paste0("/NUTSQ/NUTS3x_", year, ".csv")
  write.csv(NUTQsp,file=paste0(hydroDir,file_name))
  
  #NUTQspT=rbind(NUTQspT,NUTQsp)
  
}


stopCluster(cl)

























year=1982

for (year in yseq)
{
  # Open the NetCDF file
  print(year)
  nc <- nc_open(paste0("D:/tilloal/Documents/06_Floodrivers/HERA/dis.HERA_",year,".nc"))
  
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

# Print the first day
print(first_day)




#extract one day
#conver to raster
#aggregate by NUTS region
#repeat




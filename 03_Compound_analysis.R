
source("D:/tilloal/Documents/01_Projects/RiskDynamics/CFRisk_codes/functions_CFRisks.R")
#Set data directory
hydroDir<-("D:/tilloal/Documents/01_Projects/RiskDynamics/CFRisks_data/")


#Load primary inputs -----
outf=c()
for( Nsq in 1:88){
  print(Nsq)
  rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
  rspace=rspace[,-1]
  nrspace=rspace[Nsq,]
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

NUTS3 <- read_sf(dsn = paste0(hydroDir,"/NUTS3/NUTS3_modified.shp"))
GridNUTS3=raster( paste0(hydroDir,"/NUTS3/NUTS3_Raster3ID.tif"))
GN3=as.data.frame(GridNUTS3,xy=T)
GN3=GN3[which(!is.na(GN3[,3])),]
GN3$llcoord=paste(round(GN3$x,4),round(GN3$y,4),sep=" ")
GN3_riv=right_join(GN3,outf,by= c("llcoord"="latlong"))

GridNUTS2=raster( paste0(hydroDir,"/NUTS3/NUTS3_Raster2ID.tif"))
GN2=as.data.frame(GridNUTS2,xy=T)
GN2=GN2[which(!is.na(GN2[,3])),]
GN2$llcoord=paste(round(GN2$x,4),round(GN2$y,4),sep=" ")
GN2_riv=right_join(GN2,outf,by= c("llcoord"="latlong"))

GNF=right_join(GN3,GN2_riv,by="llcoord")

GNUTS3sf=fortify(NUTS3)

GNFx=GNF[which(is.na(GNF$NUTS3_Raster3ID)),]

#Load old and new NUTS3 region IDs
NUTS3_2010 <- read_sf(dsn = paste0(hydroDir,"/NUTS3/Regions_v2010_simplified.shp"))
NUTS3_2021 <- read_sf(dsn = paste0(hydroDir,"/NUTS3/Regions_v2021_simplified.shp"))

cores=nuts_changes[which(nuts_changes$typology=="nuts_level_3"),]


####load HANZE event set to do a first matching####
Hanze_flood=read.csv(file=paste0(hydroDir,"HANZE/HANZE_events.csv"))

#I need to create a vector of affected regions for each event
list_of_word_vectors <- lapply(Hanze_flood$Regions.affected..v2021., function(x) unlist(strsplit(x, ";")))
vectors_2021=unlist(list_of_word_vectors)
vectors_2021=unique(vectors_2021)
list_of_word_vectors2 <- lapply(Hanze_flood$Regions.affected..v2010., function(x) unlist(strsplit(x, ";")))
vectors_2010=unlist(list_of_word_vectors2)
vectors_2010=unique(vectors_2010)
N1021=cbind(vectors_2010,vectors_2021)


#load UpArea
outletname="upArea_European_01min.nc"
dir=hydroDir
outf$idlalo=paste(outf$idlo, outf$idla, sep=" ")
UpArea=UpAopen(hydroDir,outletname,outf)
head(UpArea)


#Plot parameters
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
tsize=14
osize=12
Impdates=seq(1950,2020,by=10)
valuenames=paste0("Y",Impdates)
catmap=cst7
basemap=w2

NutVector=NUTS3$NUTS_ID

#Initialize time vector
NRtime=seq.Date(as.Date("1980-01-01"),as.Date("2020-12-31"), by="day")
NRtime=as.POSIXct(NRtime)

#Load compounding hazards data
fl_events=read.csv(file=paste0(hydroDir,"/hazard_data/fl_nuts3_v1.csv"))
fl_events$enddate=as.POSIXct(fl_events$enddate)
fl_events$stardate=fl_events$enddate-(fl_events$duration*3600*24)
dr_events=read.csv(file=paste0(hydroDir,"/hazard_data/dr_nuts3_v1.csv"))
hw_events=read.csv(file=paste0(hydroDir,"/hazard_data/hw_nuts3_v1.csv"))
cw_events=read.csv(ffile=paste0(hydroDir,"/hazard_data/cw_nuts3_v1.csv"))
ws_events=read.csv(file=paste0(hydroDir,"/hazard_data/wgust_nuts3_995.csv"))
ws_events=ws_events[-which(ws_events$duration<1),]
ws_events=ws_events[,-5]


#Matching of flood events with associated hazards
Hanze_flood_1980=Hanze_flood[which(Hanze_flood$Year>1980),]
Hanze_flood_1980=Hanze_flood_1980[-which(Hanze_flood_1980$Type=="Coastal"),]
# Part 1: Detection of optimal window for flood detection with HERA high flows ----

## Initial loop exploring lag from 1 to 15 days prior to flood event ----
initF=F
if (initF==T){
  HitF=c()
  list_POTFl_large<-list()
  for (w in 1:15){
    print(w)
    wflood=c(w)
    df_POTFlood_large=c()
    Hanze_flood_M1=c()
    for (id in 1:length(NutVector)){
      print(id)
      NUTl=NutVector[id]
      #myloc=NR_Qd_total[,id]
      High_events=fl_events[which(fl_events$NUT==NUTl),]
      if (length(High_events$X)>0){
        #ms=data.frame(time=NRtime,data=myloc)
        High_events=fl_events[which(fl_events$NUT==NUTl),]
        matches <- lapply(list_of_word_vectors, function(vector) NUTl %in% vector)
        tm=which(matches==T)
        Hanze_Nut=Hanze_flood[tm,]
        Hanze_Nut=Hanze_Nut[which(Hanze_Nut$Year>1980),]
        if(length(Hanze_Nut$Country.code)>0){
          Hanze_bNut=data.frame(NUTl,Hanze_Nut)
          Hanze_flood_M1=rbind(Hanze_flood_M1,Hanze_bNut)
        }
        #now check overlap with mt detected anomalies
        df_POTFlood_keep=c()
        for (e in 1:length(Hanze_Nut$ï..ID)){
  
          #difference between end date of sim and start date of obs
          d1=(as.Date(High_events$enddate)-as.Date(Hanze_Nut$Start.date[e]))+1
          #difference between end date of sim and end date of obs
          d2=(as.Date(High_events$enddate)-as.Date(Hanze_Nut$End.date[e]))+1
          #difference between start date of sim and end date of obs
          d3=(as.Date(Hanze_Nut$End.date[e])-as.Date(High_events$stardate))-1
          #difference between start date of sim and start date of obs
          d4=(as.Date(Hanze_Nut$Start.date[e])-as.Date(High_events$stardate))-1
          
          #modelled event has to end after observed event has started
          #modelled event has to end at most xx days after oberved event has ended
          keep1=which(d1>=0 & d2<wflood[1])
          #modelled event has to start before observed event has ended
          #modelled event has to start at most xx days before observed event has started
          keep2=which(d3>=0 & d4<wflood[1])
          keep=unique(c(keep1,keep2))
          if (length(keep)>0){
            POTh_keep=High_events[keep,]
            POTh_keep$dtime=d1[keep]
            POTh_keep$dtime2=d2[keep]
            POTh_keep$dtime3=d3[keep]
            POTh_keep$dtime4=d4[keep]
            POTh_keep$NUTID=NUTl
            POTh_keep$eventID=Hanze_Nut$ï..ID[e]
            POTh_keep$eventStart=Hanze_Nut$Start.date[e]
            POTh_keep$eventEnd=Hanze_Nut$End.date[e]
            df_POTFlood_keep=rbind(df_POTFlood_keep,POTh_keep)
          } 
        }
        df_POTFlood_keep=as.data.frame(df_POTFlood_keep)
        df_POTFlood_large=rbind(df_POTFlood_large,df_POTFlood_keep)
        
      }
    }
    Hanze_flood_M1$NEID=paste0(Hanze_flood_M1$ï..ID,Hanze_flood_M1$NUTl)
    list_POTFl_large<-c(list_POTFl_large,list(df_POTFlood_large))
    
    
    #Summary of matched flood events
    Hanze_Dflood=aggregate(list(time=df_POTFlood_large$threshold),
                           by = list(NID=df_POTFlood_large$eventID),
                           FUN = function(x) c(l=length(x)))
    
    Hanze_DFmatch=inner_join(Hanze_Dflood,Hanze_flood_1980,by=c("NID"="ï..ID"))
    Hit_tot=(length(Hanze_DFmatch$NID)/length(Hanze_flood_1980$Country.code)*100)
    Hit_riv=(length(Hanze_DFmatch$NID[which(Hanze_DFmatch$Type=="River")])/length(Hanze_flood_1980$Country.code[which(Hanze_flood_1980$Type=="River")])*100)
    Hit_flash=(length(Hanze_DFmatch$NID[which(Hanze_DFmatch$Type=="Flash")])/length(Hanze_flood_1980$Country.code[which(Hanze_flood_1980$Type=="Flash")])*100)
    Hits=c(w,Hit_tot,Hit_riv,Hit_flash)
    HitF=rbind(HitF,Hits)
    
  }
  x_values=HitF[,1]
  y_values=HitF[,2]
  
  optim=elbow_finder(x_values, y_values)
  plot(HitF[,2],type="o")
  abline(v=optim[1])
}

#Optimal time lag is 4 days

## Matching of flood event with high flows with optimal lag: 4 days ----
w=4
print(w)
wflood=c(w)
df_POTFlood_large=c()
Hanze_flood_M1=c()
for (id in 1:length(NutVector)){
  print(id)
  NUTl=NutVector[id]
  myloc=NR_Qd_total[,id]
  if (length(which(!is.na(myloc)))>0){
    ms=data.frame(time=NRtime,data=myloc)
    High_events=fl_events[which(fl_events$NUT==NUTl),]
    matches <- lapply(list_of_word_vectors, function(vector) NUTl %in% vector)
    tm=which(matches==T)
    Hanze_Nut=Hanze_flood[tm,]
    Hanze_Nut=Hanze_Nut[which(Hanze_Nut$Year>1980),]
    if(length(Hanze_Nut$Country.code)>0){
      Hanze_bNut=data.frame(NUTl,Hanze_Nut)
      Hanze_flood_M1=rbind(Hanze_flood_M1,Hanze_bNut)
    }
    #now check overlap with mt detected anomalies
    df_POTFlood_keep=c()
    for (e in 1:length(Hanze_Nut$ï..ID)){
      
      #difference between end date of sim and start date of obs
      d1=(as.Date(High_events$enddate)-as.Date(Hanze_Nut$Start.date[e]))+1
      #difference between end date of sim and end date of obs
      d2=(as.Date(High_events$enddate)-as.Date(Hanze_Nut$End.date[e]))+1
      #difference between start date of sim and end date of obs
      d3=(as.Date(Hanze_Nut$End.date[e])-as.Date(High_events$stardate))-1
      #difference between start date of sim and start date of obs
      d4=(as.Date(Hanze_Nut$Start.date[e])-as.Date(High_events$stardate))-1
      
      #modelled event has to end after observed event has started
      #modelled event has to end at most xx days after oberved event has ended
      keep1=which(d1>=0 & d2<wflood[1])
      #modelled event has to start before observed event has ended
      #modelled event has to start at most xx days before observed event has started
      keep2=which(d3>=0 & d4<wflood[1])
      keep=unique(c(keep1,keep2))
      if (length(keep)>0){
        POTh_keep=High_events[keep,]
        POTh_keep$dtime=d1[keep]
        POTh_keep$dtime2=d2[keep]
        POTh_keep$dtime3=d3[keep]
        POTh_keep$dtime4=d4[keep]
        POTh_keep$NUTID=NUTl
        POTh_keep$eventID=Hanze_Nut$ï..ID[e]
        POTh_keep$eventStart=Hanze_Nut$Start.date[e]
        POTh_keep$eventEnd=Hanze_Nut$End.date[e]
        df_POTFlood_keep=rbind(df_POTFlood_keep,POTh_keep)
      } 
    }
    df_POTFlood_keep=as.data.frame(df_POTFlood_keep)
    df_POTFlood_large=rbind(df_POTFlood_large,df_POTFlood_keep)
    
  }
}
Hanze_flood_M1$NEID=paste0(Hanze_flood_M1$ï..ID,Hanze_flood_M1$NUTl)
list_POTFl_large<-c(list_POTFl_large,list(df_POTFlood_large))


#Summary of matched flood events
Hanze_Dflood=aggregate(list(time=df_POTFlood_large$threshold),
                       by = list(NID=df_POTFlood_large$eventID),
                       FUN = function(x) c(l=length(x)))

Hanze_DFmatch=inner_join(Hanze_Dflood,Hanze_flood_1980,by=c("NID"="ï..ID"))
Hit_tot=(length(Hanze_DFmatch$NID)/length(Hanze_flood_1980$Country.code)*100)
Hit_riv=(length(Hanze_DFmatch$NID[which(Hanze_DFmatch$Type=="River")])/length(Hanze_flood_1980$Country.code[which(Hanze_flood_1980$Type=="River")])*100)
Hit_flash=(length(Hanze_DFmatch$NID[which(Hanze_DFmatch$Type=="Flash")])/length(Hanze_flood_1980$Country.code[which(Hanze_flood_1980$Type=="Flash")])*100)
Hits=c(w,Hit_tot,Hit_riv,Hit_flash)

df_POTFlood=df_POTFlood_large

# Part 2: Detection of compound events with different temporal windows ----

## Definition of temporal windows ----

### Informed window----
PhCW=c(14)
PhHW=c(14)
PhWS=c(56)
PhFlood=c(56)
PhDrought=c(56)

### Three tested windows----
short=c(14)
medium=c(28)
long=c(56)


windows= c("physical", "short","medium","long")
windows= c("physical")

### Main loop over different windows ----
list_POTn_large<-list()
list_POTh_large<-list()
list_POTx_large<-list()
list_POTl_large<-list()
list_POTw_large<-list()
for (w in windows){
  print(w)
  if (w=="physical"){
    wflood=PhFlood
    wdrought=PhDrought
    wheat=PhHW
    wcold=PhCW
    wwind=PhWS
  } else if (w=="long"){
    wflood=long
    wdrought=long
    wheat=long
    wcold=long
    wwind=long
  } else if (w=="medium"){
    wflood=medium
    wdrought=medium
    wheat=medium
    wcold=medium
    wwind=medium
  } else if (w=="short"){
    wflood=short
    wdrought=short
    wheat=short
    wcold=short
    wwind=short
  }
  
  df_POTn_large=c()
  df_POTx_large=c()
  df_POTh_large=c()
  df_POTl_large=c()
  df_POTw_large=c()
  Hanze_flood_BN=c()
  for (id in 1:length(NutVector)){
    print(id)
    NUTl=NutVector[id]
    myloc=NR_Qd_total[,id]
    if (length(which(!is.na(myloc)))>0){
      ms=data.frame(time=NRtime,data=myloc)
      Low_events=dr_events[which(dr_events$NUT==NUTl),]
      High_events=fl_events[which(fl_events$NUT==NUTl),]
      HW_events=hw_events[which(hw_events$nuts_id==NUTl),]
      CW_events=cw_events[which(cw_events$nuts_id==NUTl),]
      WS_events=ws_events[which(ws_events$nuts_id==NUTl),]
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
      df_POTn_keep=c()
      df_POTw_keep=c()
      for (e in 1:length(Hanze_Nut$ï..ID)){
        #difference between end date of sim and start date of obs
        d1=(as.Date(High_events$enddate)-as.Date(Hanze_Nut$Start.date[e]))+1
        #modelled event has to end at minimum 4 days before observed starts
        #modelled event has to start at most xx days before observed event has started
        #from flood matching section
        keep1=which(d1<=(-4) & d1>(-wflood))
        keep2=c()
        keep=unique(c(keep1,keep2))
        if (length(keep)>0){
          POTh_keep=High_events[keep,]
          POTh_keep$dtime=d1[keep]
          POTh_keep$NUTID=NUTl
          POTh_keep$eventID=Hanze_Nut$ï..ID[e]
          POTh_keep$eventStart=Hanze_Nut$Start.date[e]
          POTh_keep$eventEnd=Hanze_Nut$End.date[e]
          df_POTh_keep=rbind(df_POTh_keep,POTh_keep)
        } else{
          
        }
        
        ### Drought-flood -----
        #difference between end date of sim and start date of obs
        d1=(as.Date(Low_events$endpeaks)-as.Date(Hanze_Nut$Start.date[e]))+1
        d2=(as.Date(Low_events$stpeaks)-as.Date(Hanze_Nut$End.date[e]))+1
        #modelled event has to start before observed ends
        #modelled event has to end at most xx days before observed event has started
        #from flood matching section
        keep1=which(d2<=0 & d1>(-wdrought))
        keep2=c()
        keep=unique(c(keep1,keep2))
        if (length(keep)>0){
          POTl_keep=Low_events[keep,]
          #dtime is the difference between end date of pre and start date of flood
          POTl_keep$dtime=d1[keep]
          POTl_keep$NUTID=NUTl
          POTl_keep$eventStart=Hanze_Nut$Start.date[e]
          POTl_keep$eventEnd=Hanze_Nut$End.date[e]
          POTl_keep$eventID=Hanze_Nut$ï..ID[e]
          df_POTl_keep=rbind(df_POTl_keep,POTl_keep)
        }
        
        ### Hot-Wet sequence -----
        
        #difference between end date of sim and start date of obs
        d1=(as.Date(HW_events$end)-as.Date(Hanze_Nut$Start.date[e]))+1
        # d2=(as.Date(HW_events$end)-as.Date(Hanze_Nut$End.date[e]))+1
        #modelled event has to end before observed starts
        #modelled event has to end at most xx days before observed event starts
        #from flood matching section
        keep1=which(d1<=0 & d1>(-wheat))
        keep2=c()
        keep=unique(c(keep1,keep2))
        if (length(keep)>0){
          POTx_keep=HW_events[keep,]
          POTx_keep$dtime=d1[keep]
          POTx_keep$NUTID=NUTl
          POTx_keep$eventStart=Hanze_Nut$Start.date[e]
          POTx_keep$eventEnd=Hanze_Nut$End.date[e]
          POTx_keep$eventID=Hanze_Nut$ï..ID[e]
          df_POTx_keep=rbind(df_POTx_keep,POTx_keep)
        }
        
        ### Cold-flood event -----
        #difference between end date of sim and start date of obs
        d1=(as.Date(CW_events$end)-as.Date(Hanze_Nut$Start.date[e]))+1
        d2=(as.Date(CW_events$begin)-as.Date(Hanze_Nut$End.date[e]))+1
        #modelled event has to start before observed ends
        #modelled event has to end at most xx days before observed event starts
        #from flood matching section
        keep1=which(d2<=0 & d1>(-wcold))
        keep2=c()
        keep=unique(c(keep1,keep2))
        if (length(keep)>0){
          POTn_keep=CW_events[keep,]
          POTn_keep$dtime=d1[keep]
          POTn_keep$NUTID=NUTl
          POTn_keep$eventID=Hanze_Nut$ï..ID[e]
          POTn_keep$eventEnd=Hanze_Nut$End.date[e]
          POTn_keep$eventStart=Hanze_Nut$Start.date[e]
          df_POTn_keep=rbind(df_POTn_keep,POTn_keep)
        }

        ### Wind-flood event -----
        #difference between end date of sim and start date of obs
        d1=(as.Date(WS_events$end)-as.Date(Hanze_Nut$Start.date[e]))+1
        d2=(as.Date(WS_events$begin)-as.Date(Hanze_Nut$End.date[e]))+1
        #modelled event has to start before observed ends
        #modelled event has to end at most xx days before observed event starts
        #from flood matching section
        keep1=which(d2<=0 & d1>(-wwind))
        keep2=c()
        keep=unique(c(keep1,keep2))
        if (length(keep)>0){
          POTw_keep=WS_events[keep,]
          POTw_keep$dtime=d1[keep]
          POTw_keep$NUTID=NUTl
          POTw_keep$eventID=Hanze_Nut$ï..ID[e]
          POTw_keep$eventEnd=Hanze_Nut$End.date[e]
          POTw_keep$eventStart=Hanze_Nut$Start.date[e]
          df_POTw_keep=rbind(df_POTw_keep,POTw_keep)
        }
      }
      
      df_POTh_keep=as.data.frame(df_POTh_keep)
      df_POTl_keep=as.data.frame(df_POTl_keep)
      df_POTx_keep=as.data.frame(df_POTx_keep)
      df_POTn_keep=as.data.frame(df_POTn_keep)
      df_POTw_keep=as.data.frame(df_POTw_keep)
      
      df_POTh_large=rbind(df_POTh_large,df_POTh_keep)
      df_POTl_large=rbind(df_POTl_large,df_POTl_keep)
      df_POTx_large=rbind(df_POTx_large,df_POTx_keep)
      df_POTn_large=rbind(df_POTn_large,df_POTn_keep)
      df_POTw_large=rbind(df_POTw_large,df_POTw_keep)
      
    }
  }
  
  Hanze_flood_BN$NEID=paste0(Hanze_flood_BN$ï..ID,Hanze_flood_BN$NUTl)
  df_POTh_large$NEID=paste0(df_POTh_large$eventID,df_POTh_large$NUTID)
  df_POTl_large$NEID=paste0(df_POTl_large$eventID,df_POTl_large$NUTID)
  df_POTx_large$NEID=paste0(df_POTx_large$eventID,df_POTx_large$NUTID)
  df_POTn_large$NEID=paste0(df_POTn_large$eventID,df_POTn_large$NUTID)
  df_POTw_large$NEID=paste0(df_POTw_large$eventID,df_POTw_large$NUTID)
  
  list_POTn_large<-c(list_POTn_large,list(df_POTn_large))
  list_POTx_large<-c(list_POTx_large,list(df_POTx_large))
  list_POTh_large<-c(list_POTh_large,list(df_POTh_large))
  list_POTl_large<-c(list_POTl_large,list(df_POTl_large))
  list_POTw_large<-c(list_POTw_large,list(df_POTw_large))
  
}

## Analysis for multiple windows ----
plot1=list()
list_eventc2=list()
list_call=list()
for (w in c(1:4)){
  print(w)
  Events_prefflood=inner_join(list_POTh_large[[w]],Hanze_flood_1980, by=c("eventID"="ï..ID"))
  
  #verification: preflood with flood events
  df_POTFlood$NEID=paste0(df_POTFlood$NUTID,df_POTFlood$eventID)
  matF=na.omit(match(df_POTFlood$NEID,Events_prefflood$NEID))
  
  #second filter for prefloods based on start date
  dX=difftime(as.Date(Events_prefflood$stardate),as.Date(Events_prefflood$Start.date),units="day")
  Events_prefflood$dtime2=dX
  
  #aggregate by region to have a kind of success rate spatially
  Hanze_SDflood=aggregate(list(he=df_POTFlood$threshold),
                          by = list(NUT=df_POTFlood$NUT,eid=df_POTFlood$eventID),
                          FUN = function(x) c(l=length(x)))
  
  Hanze_SDflood=aggregate(list(he=Hanze_SDflood$he),
                          by = list(NUT=Hanze_SDflood$NUT),
                          FUN = function(x) c(l=length(x)))
  
  Hanze_SOflood=aggregate(list(he=Hanze_flood_BN$Year),
                          by = list(NUT=Hanze_flood_BN$NUTl),
                          FUN = function(x) c(l=length(x)))
  
  Hanze_spRat=full_join(Hanze_SOflood,Hanze_SDflood,by="NUT")
  Hanze_spRat$succerat=Hanze_spRat$he.y/Hanze_spRat$he.x*100
  
  Hanze_spRat=inner_join(Hanze_spRat,GNUTS3sf, by=c("NUT"="NUTS_ID"))
 
  
  legend="Detection rate (%)"
  palet=c(hcl.colors(11, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
  ggplot(basemap) +
    geom_sf(fill="white",color="darkgrey",size=0.5)+
    geom_sf(data=Hanze_spRat,aes(fill=succerat,geometry=geometry),alpha=0.9,color="transparent")+
    # scale_fill_manual(
    #   values=palet, breaks=rev(labels), name=legend2)   +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_fill_gradientn(
      colors=palet,
      oob = scales::squish,na.value="transparent", name=legend)   +
    labs(x="Longitude", y = "Latitude")+
    #guides(fill = guide_coloursteps(barwidth = 1.5, barheight = 14,nbins=5))+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=osize),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))
  
  ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/MH_events_8120_HANZE.jpg"), width=30, height=20, units=c("cm"),dpi=300) 
  
  #aggregate by event ID
  Hanze_flood_BN$NEID=paste0(Hanze_flood_BN$ï..ID,Hanze_flood_BN$NUTl)
  
  #Matched flood events
  Ev_flmatch=aggregate(list(time=df_POTFlood$time),
                       by = list(NID=df_POTFlood$NUTID,EID=df_POTFlood$eventID),
                       FUN = function(x) c(l=length(x)))
  names(Ev_flmatch)[3]="N_match"
  Ev_flmatch$NEID=paste0(Ev_flmatch$EID,Ev_flmatch$NID)
  Ev_flood1=full_join(Ev_flmatch,Hanze_flood_BN,by="NEID")
  
  
  #Precondition flood events
  Ev_flmatch=aggregate(list(time=Events_prefflood$time),
                       by = list(NID=Events_prefflood$NUTID,EID=Events_prefflood$eventID),
                       FUN = function(x) c(l=length(x)))
  names(Ev_flmatch)[3]="N_PreEv"
  Ev_flmatch$NEID=paste0(Ev_flmatch$EID,Ev_flmatch$NID)
  Ev_flood2=full_join(Ev_flmatch,Hanze_flood_BN,by="NEID")
  
  
  Events_precoflood=inner_join(list_POTl_large[[w]],Hanze_flood_1980, by=c("eventID"="ï..ID"))
  Events_prehflood=inner_join(list_POTx_large[[w]],Hanze_flood_1980, by=c("eventID"="ï..ID"))
  Events_precflood=inner_join(list_POTn_large[[w]],Hanze_flood_1980, by=c("eventID"="ï..ID"))
  Events_prewflood=inner_join(list_POTw_large[[w]],Hanze_flood_1980, by=c("eventID"="ï..ID"))
  
  Events_prehflood$time=Events_prehflood$begin
  Events_precflood$time=Events_precflood$begin
  Events_prewflood$time=Events_prewflood$begin
  
  #creation of table with similar structures for each hazard
  C_flood_ev=Events_prefflood
  C_flood_ev$stdate=as.POSIXct(as.Date(C_flood_ev$stpeaks, origin = "1980-01-01"))
  colnames(C_flood_ev)
  #identification of events that correspond to the recorded flood event
  #needed colums:  "NUT"      "dtime"    "duration" "peaks"    "time"     "eventID"  "NEID"  
  C_flood_ev=C_flood_ev[c(10,12,8,3,7,14,17)]
  C_flood_ev$peaks=(C_flood_ev$peaks-.98)/0.02
  hist(C_flood_ev$peaks)
  C_flood_ev$time=as.Date(C_flood_ev$time)

  C_cold_ev=Events_precflood
  dura=as.numeric(difftime(as.Date(C_cold_ev$end),as.Date(C_cold_ev$begin)))
  C_cold_ev$time=as.Date(C_cold_ev$time)
  C_cold_ev$duration=C_cold_ev$intensity
  C_cold_ev$NUT=C_cold_ev$nuts_id
  C_cold_ev$peaks=1
  C_cold_ev$dtime[which(C_cold_ev$dtime>0)]=0
  colsel=match(colnames(C_flood_ev),colnames(C_cold_ev))
  C_cold_ev=C_cold_ev[colsel]
  ymon=month(C_cold_ev$time)
  #remove cold events occuring during the "warm season"
  season=rep("cold",length(ymon))
  season[which(ymon>4 & ymon<11)]="warm"
  C_cold_ev=C_cold_ev[-which(season=="warm"),]
  
  C_heat_ev=Events_prehflood
  dura=as.numeric(difftime(as.Date(C_heat_ev$end),as.Date(C_heat_ev$begin)))
  C_heat_ev$time=as.Date(C_heat_ev$time)
  C_heat_ev$duration=C_heat_ev$intensity
  C_heat_ev$NUT=C_heat_ev$nuts_id
  C_heat_ev$peaks=1
  colsel=match(colnames(C_flood_ev),colnames(C_heat_ev))
  C_heat_ev=C_heat_ev[colsel]
  ymon=month(C_heat_ev$time)
  #remove hot events occuring during the "cold season"
  season=rep("cold",length(ymon))
  season[which(ymon>4 & ymon<11)]="warm"
  C_heat_ev=C_heat_ev[-which(season=="cold"),]
  
  
  C_wind_ev=Events_prewflood
  dura=as.numeric(difftime(as.Date(C_wind_ev$end),as.Date(C_wind_ev$begin)))
  C_wind_ev$time=as.Date(C_wind_ev$time)
  C_wind_ev$duration=dura
  C_wind_ev$NUT=C_wind_ev$nuts_id
  C_wind_ev$peaks=C_wind_ev$intensity
  C_wind_ev$dtime[which(C_wind_ev$dtime>0)]=0
  colsel=match(colnames(C_flood_ev),colnames(C_wind_ev))
  C_wind_ev=C_wind_ev[colsel]
  
  C_drought_ev=Events_precoflood
  C_drought_ev$low_peaks=((1-C_drought_ev$low_peaks)-0.95)/0.05
  C_drought_ev$time=as.Date(C_drought_ev$time)
  C_drought_ev$peaks=C_drought_ev$low_peaks
  C_drought_ev$duration=C_drought_ev$cluster_length
  C_drought_ev$dtime[which(C_drought_ev$dtime>0)]=0
  colsel=match(colnames(C_flood_ev),colnames(C_drought_ev))
  C_drought_ev=C_drought_ev[colsel]
  names(C_drought_ev)=names(C_heat_ev)=names(C_cold_ev)=names(C_wind_ev)=names(C_flood_ev)
  
  C_drought_ev$evtype="drought"
  C_flood_ev$evtype="flood"
  C_heat_ev$evtype="heat"
  C_cold_ev$evtype="cold"
  C_wind_ev$evtype="wind"
  
  C_all_haz=rbind(C_flood_ev,C_drought_ev,C_heat_ev,C_cold_ev,C_wind_ev)
  list_call<-c(list_call,list(C_all_haz))
  
  all_events=unique(Hanze_flood_1980$ï..ID)
  mh_events=unique(C_all_haz$eventID)
  
  #number of events for each category
  v1=(unique(Hanze_flood_1980$ï..ID))
  v2=(unique(Events_prefflood$eventID))
  v3=(unique(Events_precoflood$eventID))
  v4=(unique(Events_prehflood$eventID))
  v5=(unique(Events_precflood$eventID))
  v6=(unique(Events_prewflood$eventID))
  
  #extract the number of single events
  s_events=all_events[which(is.na(match(all_events,mh_events)))]
  
  length(all_events)
  length(mh_events)
  length(s_events)
  
  evtypes=c("drought","flood","heat","cold","wind","single")
  events_class=c()
  for (ie in 1:length(all_events))
  {
    evsel=all_events[ie]
    cls=vector(mode="integer",length=7)
    C_all_s=C_all_haz[which(C_all_haz$eventID==evsel),]
    cls[1]=evsel
    evt=unique(C_all_s$evtype)
    if (length(C_all_s$NUT)<1){
      evt="single"
    }
    mt=which(!is.na(match(evtypes,evt)))
    cls[mt+1]=1
    events_class=rbind(events_class,cls)
  }
  
  events_class=as.data.frame(events_class)
  events_class$nh=events_class[,2]+events_class[,3]+events_class[,4]+events_class[,5]+events_class[,6]
  
  
  colnames(events_class)=c("eventID","predrought","preflood","preheat","precold","prewind","single","nh")
  
  #event categories
  mh_cat=c("single","drought-flood","wet-sequence","heat-flood","cold-flood","heat-drought-flood","cold-wet-sequence","wind-flood")
  
  #Classification into all possible combinations of hazards
  
  #single events category
  sing_events=events_class[which(events_class$nh==0),]
  sing_events$event_name="single-flood"
  
  #bivariate events category
  biv_events=events_class[which(events_class$nh==1),]
  biv_events$event_name=NA
  biv_events$event_name[which(biv_events$predrought==1)]="drought-flood"
  biv_events$event_name[which(biv_events$preflood==1)]="wet-sequence"
  biv_events$event_name[which(biv_events$preheat==1)]="heat-flood"
  biv_events$event_name[which(biv_events$precold==1)]="cold-flood"
  biv_events$event_name[which(biv_events$prewind==1)]="wind-flood"
  
  #trvariate events category
  triv_events=events_class[which(events_class$nh==2),]
  triv_events$event_name=NA
  triv_events$event_name[which(triv_events$predrought==1 & triv_events$preflood==1)]="drought-wet-sequence"
  triv_events$event_name[which(triv_events$predrought==1 & triv_events$preheat==1)]="heat-drought-flood"
  triv_events$event_name[which(triv_events$preflood==1 & triv_events$preheat==1)]="heat-wet-sequence"
  triv_events$event_name[which(triv_events$preflood==1 & triv_events$precold==1)]="cold-wet-sequence"
  triv_events$event_name[which(triv_events$preheat==1 & triv_events$precold==1)]="cold-heat-wet"
  triv_events$event_name[which(triv_events$predrought==1 & triv_events$precold==1)]="cold-drought-flood"
  triv_events$event_name[which(triv_events$predrought==1 & triv_events$prewind==1)]="wind-drought-flood"
  triv_events$event_name[which(triv_events$preflood==1 & triv_events$prewind==1)]="wind-wet-sequence"
  triv_events$event_name[which(triv_events$precold==1 & triv_events$prewind==1)]="wind-cold-flood"
  triv_events$event_name[which(triv_events$preheat==1 & triv_events$prewind==1)]="wind-heat-flood"
  
  
  #all events together
  event_class2=rbind(sing_events,biv_events,triv_events)
  
  #quadrivariate events category
  quadriv_events=events_class[which(events_class$nh==3),]
  if (length(quadriv_events$eventID)>0){
    quadriv_events$event_name=NA
    quadriv_events$event_name[which(quadriv_events$predrought==1 & quadriv_events$preflood==1 & quadriv_events$preheat==1)]="heat-drought-to-wet-sequence"
    quadriv_events$event_name[which(quadriv_events$predrought==1 & quadriv_events$preflood==1 & quadriv_events$precold==1)]="cold-drought-to-wet-sequence"
    quadriv_events$event_name[which(quadriv_events$preheat==1 & quadriv_events$preflood==1 & quadriv_events$precold==1)]="cold-heat-wet-sequence"
    quadriv_events$event_name[which(quadriv_events$prewind==1 & quadriv_events$preheat==1 & quadriv_events$precold==1)]="wind-cold-heat-flood"
    quadriv_events$event_name[which(quadriv_events$predrought==1 & quadriv_events$prewind==1 & quadriv_events$precold==1)]="wind-cold-drought-flood"
    quadriv_events$event_name[which(quadriv_events$predrought==1 & quadriv_events$preheat==1 & quadriv_events$prewind==1)]="wind-heat-drought-flood"
    quadriv_events$event_name[which(quadriv_events$preflood==1 & quadriv_events$preheat==1 & quadriv_events$prewind==1)]="wind-heat-wet-sequence"
    quadriv_events$event_name[which(quadriv_events$preflood==1 & quadriv_events$precold==1 & quadriv_events$prewind==1)]="wind-cold-wet-sequence"
    quadriv_events$event_name[which(quadriv_events$preflood==1 & quadriv_events$predrought==1 & quadriv_events$prewind==1)]="wind-drought-to-wet-sequence"
    quadriv_events$event_name[which(quadriv_events$precold==1 & quadriv_events$preheat==1 & quadriv_events$predrought==1)]="cold-heat-drought-flood"
    event_class2=rbind(event_class2,quadriv_events)
  }
  #remaining events
  quindriv_events=events_class[which(events_class$nh==4),]
  if (length(quindriv_events$eventID)>0){
    quindriv_events$event_name=NA
    quindriv_events$event_name[which(quindriv_events$precold==1,quindriv_events$predrought==1 &
                                       quindriv_events$preflood==1 & quindriv_events$preheat==1)]="heat-cold-drought-to-wet-sequence"
    quindriv_events$event_name[which(quindriv_events$precold==1,quindriv_events$predrought==1 &
                                       quindriv_events$preflood==1 & quindriv_events$prewind==1)]="wind-cold-drought-to-wet-sequence"
    quindriv_events$event_name[which(quindriv_events$prewind==1,quindriv_events$predrought==1 & 
                                       quindriv_events$preflood==1 & quindriv_events$preheat==1)]="wind-heat-drought-to-wet-sequence"
    quindriv_events$event_name[which(quindriv_events$prewind==1,quindriv_events$precold==1 & 
                                       quindriv_events$preflood==1 & quindriv_events$preheat==1)]="wind-heat-cold-wet-sequence"
    quindriv_events$event_name[which(quindriv_events$prewind==1,quindriv_events$predrought==1 & 
                                       quindriv_events$precold==1 & quindriv_events$preheat==1)]="wind-heat-cold-drought-flood"
    event_class2=rbind(event_class2,quindriv_events)
  }
  
  sixdriv_events=events_class[which(events_class$nh==5),]
  if (length(sixdriv_events$eventID)>0){
    sixdriv_events$event_name=NA
    sixdriv_events$event_name="wind-heat-cold-drought-to-wet-sequence"
    
    event_class2=rbind(event_class2,sixdriv_events)
  }
  
  uni=unique(event_class2$event_name)
  event_class2$event_name <- factor(event_class2$event_name, 
                                    levels =uni)
  event_class2$nhcat=as.character(event_class2$nh)
  list_eventc2=c(list_eventc2,list(event_class2))
  
  palet=c(hcl.colors(5, palette = "YlOrRd", alpha = NULL, rev = T, fixup = TRUE))
  ###[Supplement Plot] - plot pf detected multi-hazard events before attribution to archetypes
  
  plot1[[w]]<- ggplot(event_class2, aes(x = event_name,fill = nhcat))+ 
    geom_bar() + 
    coord_flip()+
    #scale_fill_gradientn(colors=palet,name="Number of compounding hazards") +
    scale_fill_manual(
      values=palet, name="Number of compounding hazards")   +
    # guides(fill = guide_colorsteps(barwidth = 20, barheight = 1))+
    guides(fill = guide_legend(title.position = "top",nrow = 1,
                               reverse=F, title.hjust=0.5, 
                               barwidth = 3, barheight = 1.5)) +
    scale_y_continuous(name="Number of events",expand = expansion(mult = c(0, 0)))+
    scale_x_discrete(name="Compound event name")+
    theme(axis.title=element_text(size=14, face="bold"),
          axis.text = element_text(size=12),
          axis.text.x = element_text(size=12,face="bold"),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=12),
          axis.ticks.y = element_blank(),  # Remove x-axis tickmarks
          panel.grid.major.y = element_line(color = "lightgray",linetype="dashed"),
          panel.grid.minor.y = element_line(color = "lightgray"),
          panel.grid.major.x = element_blank(),
          legend.position = "bottom",
          legend.title.position = "top",
          legend.text.position = "bottom",
          # panel.grid.major = element_line(colour = "grey80"),
          panel.grid.minor.x = element_blank(),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/compound_events_hanzeFX_",w,".jpg"),  plot1[[w]], width=30, height=20, units=c("cm"),dpi=300) 
  
}


###For several window assessment -----
newdata=list_call[[1]]
nev=unique(newdata$eventID)
#aggregate by hazard
all_events=unique(Hanze_flood_1980$ï..ID)
all_events=unique(Hanze_flood_BN$NEID)
evtypes=c("drought","flood","heat","cold","wind")
events_class=c()
for (evi in 1:length(all_events)){
  evid=all_events[evi]
  #evid="2413ES521"
  myev=newdata[which(!is.na(match(newdata$NEID,evid))),]
  uh=unique(myev$evtype)
  cls=rep(NA, 17)
  cls[1]=evid
  evt=uh
  if (length(myev$NUT)>0){
    cls[2]=Hanze_flood_BN$NUTl[evi]
  }
  for (hz in unique(myev$evtype)){
    myehz=myev[which(myev$evtype==hz),]
    #grouping by lag
    myehz$grp=1
    myehz=myehz[order(myehz$dtime),]
    grp=1
    if (length(myehz$NUT)>1){
      for (dt in 2:length(myehz$dtime)){
        dti=myehz$dtime[dt-1]-myehz$dtime[dt]
        if(dti<=(-2)){
          grp=grp+1
        }
        myehz$grp[dt]=grp
      }
    }
    lev=length(unique(myehz$grp))
    mxp=round(max(myehz$peaks),3)
    md=which.min(abs(myehz$dtime))
    if (length(myev$NUT)<1){
      evt="single"
    }
    mt=which(!is.na(match(evtypes,hz)))
    cls[(mt-1)*3+3]=lev
    cls[(mt-1)*3+4]=max(sum(myehz$peaks*myehz$duration),0)
    cls[(mt-1)*3+5]=myehz$dtime[md]
  }
  events_class=rbind(events_class,cls)
} 
events_class=as.data.frame(events_class)
colnames(events_class)=c("eventxNUT_ID","NUTSID","drought.n","drought.i","drought.lag",
                         "flood.n","flood.i","flood.lag",
                         "HW.n","HW.i","HW.lag",
                         "CW.n","CW.i","CW.lag",
                         "WS.n","WS.i","WS.lag")
events_class[,c(3:17)]=lapply( events_class[,c(3:17)], function(x) as.numeric(as.character(x)))


#match with Hanze
#remove coastal flooding events
Hanze_flood_BN2=Hanze_flood_BN[-which(Hanze_flood_BN$Type=="Coastal"),]
event_classAI=inner_join(events_class,Hanze_flood_BN2,by=c("eventxNUT_ID"="NEID"))

write.csv(event_classAI,file=paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/CompoundEventsXHanze_v3.csv"))


#Add RP of the flood event to previous dataset
Hanze_events_RP=read.csv(file="D:/tilloal/Documents/01_Projects/RiskDynamics/Data/Attribution_losses_under_factual_scenario_river_with_RP.csv")
Hanze_events_RP$NEID=paste0(Hanze_events_RP$HANZE_ID,Hanze_events_RP$NUTS3)
unikid=unique(Hanze_events_RP$HANZE_ID)
Hanze_events_RP2=c()

for (evx in unikid){
  sev=Hanze_events_RP[which(Hanze_events_RP$HANZE_ID==evx),]
  sumloss=sum(sev$Economic_loss)
  sev$ratioloss=round(sev$Economic_loss/sumloss,2)
  
  sumfat=sum(sev$Fatalities)
  sev$ratiofat=round(sev$Fatalities/sumfat,2)
  
  sumaf=sum(sev$Persons_affected)
  sev$ratioaf=round(sev$Persons_affected/sumaf,2)
  
  Hanze_events_RP2=rbind(Hanze_events_RP2,sev)
}

Hanze_events_RP2=Hanze_events_RP2[,-c(1,10,11,12)]
rpm=match(Hanze_flood_BN$NEID,Hanze_events_RP$NEID)
Hanze_flood_BN$RPflood=Hanze_events_RP$NEID[rpm]

taglist=c()
event_classic=c()
for (w in 2:4){
  event_all1=list_eventc2[[w]]
  aev=length(unique(event_all1$eventID))
  event_mh=event_all1[which(event_all1$nh>=2),]
  mhev=unique(event_mh$eventID)
  event_nh=event_all1[which(event_all1$nh<2),]
  #loop over event_mh
  event_m1=event_mh[1,]
  # Repeat each row nh times
  repeated_table <- event_mh[rep(1:nrow(event_mh), event_mh$nh), ]
  event_mh2=c()
  for (evi in 1:length(mhev)){
    event_m1=repeated_table[which(repeated_table$eventID==mhev[evi]),]
    hs=which(event_m1[1,]==1)
    locp=c(2,3,4,5,6)
    for(r in 1:length(hs)){
      rh=hs[r]
      locn=locp[-which(locp==rh)]
      event_m1[r,locn]=0
    }
    event_mh2=rbind(event_mh2,event_m1)
  }
  
  event_mh2$event_name[which(event_mh2$predrought==1)]="drought-flood"
  event_mh2$event_name[which(event_mh2$preflood==1)]="wet-sequence"
  event_mh2$event_name[which(event_mh2$preheat==1)]="heat-flood"
  event_mh2$event_name[which(event_mh2$precold==1)]="cold-flood"
  event_mh2$event_name[which(event_mh2$prewind==1)]="windstorm-flood"
  
  #remerge events
  event_mh3=rbind(event_nh,event_mh2)
  
  #merge with HANZE
  event_class3=inner_join(event_mh3,Hanze_flood_1980,by=c("eventID"="ï..ID"))
  event_classic=rbind(event_classic,event_class3)
  
  #aggregate by event type for fatalities, people affected and economic damages
  event_tag<-aggregate(list(losses=event_class3$Losses..2020.euro.,fat=event_class3$Fatalities,af=event_class3$Persons.affected),
                       by = list(Events=event_class3$event_name),
                       FUN = function(x) c(mean=mean(x,na.rm=T),len=length(x),sum=sum(x,na.rm=T)))
  event_tag <- do.call(data.frame, event_tag)
  
  
  un_ev=event_tag$Events
  matnh=match(un_ev,event_class3$event_name)
  event_tag$nh=event_class3$nh[matnh]
  event_tag$window=windows[w]
  taglist=rbind(taglist,event_tag)
}

color_events=c("single-flood"="gray","wet-sequence"="dodgerblue4","heat-flood"="darksalmon","drought-flood"="orange","cold-flood"="lightblue","windstorm-flood"="seagreen")

color_windows=c("short"="darksalmon","medium"="orange","long"="lightblue","physical"="darkgrey")

# reorder the levels of the window column
taglist$Events <- factor(taglist$Events, 
                         levels = c("single-flood","wet-sequence","heat-flood","drought-flood","cold-flood","windstorm-flood"),
                         ordered = TRUE)

taglist$window <- factor(taglist$window, 
                         levels =c("short","medium","long","physical"),
                         ordered = TRUE)

taghigh=taglist[(c(4,7,9,20,23,24)),]
taghigh$color="black"
ggplot(taglist[-c(1,2,3,4,5,6),], aes(x = Events,y=losses.len,fill=window,category=window))+ 
  geom_bar(stat = "identity",position="dodge",color="black")+ 
  geom_bar(data=taghigh, aes(x = Events,y=losses.len,color=color),
           stat = "identity",position="dodge",fill="transparent",lwd=1,lty=5)+ 
  scale_y_continuous(name="number of events")+
  coord_flip()+
  scale_fill_manual(
    values=color_windows,name="windows") +
  scale_color_manual(
    values="black",label="",name="informed \nwindow") +
  theme(axis.title=element_text(size=18, face="bold"),
        axis.text = element_text(size=14),
        axis.text.x = element_text(size=14,face="bold"),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=18, face="bold"),
        legend.text = element_text(size=14),
        axis.ticks.y = element_blank(),  # Remove x-axis tickmarks
        panel.grid.major.y = element_line(color = "lightgray"),
        panel.grid.minor.y = element_line(color = "lightgray"),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/compound_events_windows_newF.jpg"), width=30, height=20, units=c("cm"),dpi=300) 



# Part 3: Analysis on the informed window ----
w=1
##Event-based analysis ----
#create final table with compopund events archetypes
event_all1=list_eventc2[[w]]
aev=length(unique(event_all1$eventID))
event_mh=event_all1[which(event_all1$nh>=2),]
mhev=unique(event_mh$eventID)
event_nh=event_all1[which(event_all1$nh<2),]

event_m1=event_mh[1,]

# Repeat each row nh times
repeated_table <- event_mh[rep(1:nrow(event_mh), event_mh$nh), ]
event_mh2=c()
for (evi in 1:length(mhev)){
  event_m1=repeated_table[which(repeated_table$eventID==mhev[evi]),]
  hs=which(event_m1[1,]==1)
  locp=c(2,3,4,5,6)
  for(r in 1:length(hs)){
    rh=hs[r]
    locn=locp[-which(locp==rh)]
    event_m1[r,locn]=0
  }
  event_mh2=rbind(event_mh2,event_m1)
}

event_mh2$event_name[which(event_mh2$predrought==1)]="drought-flood"
event_mh2$event_name[which(event_mh2$preflood==1)]="wet-sequence"
event_mh2$event_name[which(event_mh2$preheat==1)]="heat-flood"
event_mh2$event_name[which(event_mh2$precold==1)]="cold-flood"
event_mh2$event_name[which(event_mh2$prewind==1)]="windstorm-flood"
#remerge events
event_mh3=rbind(event_nh,event_mh2)
length(which(event_mh3$nh==0))/length(event_mh3$eventID)
#merge with HANZE
event_class3=inner_join(event_mh3,Hanze_flood_1980,by=c("eventID"="ï..ID"))

#aggregate by events and years
event_tag<-aggregate(list(losses=event_class3$Losses..2020.euro.,fat=event_class3$Fatalities,af=event_class3$Persons.affected),
                     by = list(Events=event_class3$event_name,year=event_class3$Year),
                     FUN = function(x) c(mean=mean(x,na.rm=T),len=length(x),sum=sum(x,na.rm=T)))
event_tag <- do.call(data.frame, event_tag)

mylab= "number of events"


## Region-event based analysis ----
event_allR1=event_classAI
aev=length(unique(event_allR1$eventxNUT_ID))
#1 if the hazard occured, otherwise 0
event_allR1$dn=0
event_allR1$dn[which(!is.na(event_allR1$drought.n))]=1
event_allR1$fn=0
event_allR1$fn[which(!is.na(event_allR1$flood.n))]=1
event_allR1$hn=0
event_allR1$hn[which(!is.na(event_allR1$HW.n))]=1
event_allR1$cn=0
event_allR1$cn[which(!is.na(event_allR1$CW.n))]=1
event_allR1$wn=0
event_allR1$wn[which(!is.na(event_allR1$WS.n))]=1

event_allR1$nh=event_allR1$dn + event_allR1$fn + event_allR1$hn + event_allR1$cn + event_allR1$wn
event_mhR=event_allR1[which(event_allR1$nh>=2),]
mhev=unique(event_mhR$eventxNUT_ID)
event_mhR=event_mhR[,c(1,2,18:25,39:44)]
# mhev/aev
aev=length(unique(event_allR1$eventxNUT_ID))
event_nhR=event_allR1[which(event_allR1$nh<2),]
event_nhR=event_nhR[,c(1,2,18:25,39:44)]
#loop over event_mh

# Repeat each row nh times
repeated_table <- event_mhR[rep(1:nrow(event_mhR), event_mhR$nh), ]
event_mhR2=c()
for (evi in 1:length(mhev)){
  event_r1=repeated_table[which(repeated_table$eventxNUT_ID==mhev[evi]),]
  hs=which(event_r1[1,]==1)
  locp=c(11,12,13,14,15)
  for(r in 1:length(hs)){
    rh=hs[r]
    locn=locp[-which(locp==rh)]
    event_r1[r,locn]=0
  }
  event_mhR2=rbind(event_mhR2,event_r1)
}

event_mhR2$event_name=NA
event_mhR2$event_name[which(event_mhR2$dn==1)]="drought-flood"
event_mhR2$event_name[which(event_mhR2$fn==1)]="wet-sequence"
event_mhR2$event_name[which(event_mhR2$hn==1)]="heat-flood"
event_mhR2$event_name[which(event_mhR2$cn==1)]="cold-flood"
event_mhR2$event_name[which(event_mhR2$wn==1)]="windstorm-flood"
#remerge events
event_nhR$event_name="single-flood"
event_nhR$event_name[which(event_nhR$dn==1)]="drought-flood"
event_nhR$event_name[which(event_nhR$fn==1)]="wet-sequence"
event_nhR$event_name[which(event_nhR$hn==1)]="heat-flood"
event_nhR$event_name[which(event_nhR$cn==1)]="cold-flood"
event_nhR$event_name[which(event_nhR$wn==1)]="windstorm-flood"
event_mhR3=rbind(event_nhR,event_mhR2)


event_MHRn<-aggregate(list(mh=event_mhR2$nh),
                    by = list(NUT=event_mhR2$NUTl),
                    FUN = function(x) c(len=length(x)))
event_MHRn <- do.call(data.frame, event_MHRn)

#map from event_MHRn: total number of mh events per NUTS
event_Nt<-aggregate(list(nh=event_mhR3$NUTl),
                     by = list(NUT=event_mhR3$NUTl,event_type=event_mhR3$event_name),
                     FUN = function(x) c(len=length(x)))
event_Nt <- do.call(data.frame, event_Nt)

unique(event_Nt$event_type)

#map from event_Nt2: number of unique mh event per NUTS
event_Nt2<-aggregate(list(ev=event_Nt$event_type),
                    by = list(NUT=event_Nt$NUT),
                    FUN = function(x) c(len=length(x)))
event_Nt2 <- do.call(data.frame, event_Nt2)

event_Nt3<-aggregate(list(ev=event_Nt$nh),
                     by = list(NUT=event_Nt$event_type),
                     FUN = function(x) c(sum=sum(x)))
event_Nt3 <- do.call(data.frame, event_Nt3)

event_Nt2$dominant=NA
for (ii in 1:length(event_Nt2$NUT)){
  evx=event_Nt[which(event_Nt$NUT==event_Nt2$NUT[[ii]]),]
  evmax=which.max(evx$nh)
  mhx=evx$event_type[evmax]
  event_Nt2$dominant[ii]=mhx
}

length(which(event_Nt2$dominant=="single-flood"))/length(event_Nt2$dominant)
event_NtX=aggregate(list(ev=event_Nt2$NUT),
                    by = list(NUT=event_Nt2$dominant),
                    FUN = function(x) c(len=length(x)))
event_NtX <- do.call(data.frame, event_NtX)
mylab= "number of events"



### [Plot] - Evolution of the frequency of compound events through time ----
unev=unique(event_tag$Events)
#fill gaps for years with no recorded events
all_years=data.frame(year = rep(seq(min(event_tag$year), max(event_tag$year), by = 1),6))
all_evt=data.frame(event=c(rep(unev[1],40),rep(unev[2],40),rep(unev[3],40),rep(unev[4],40),rep(unev[5],40),rep(unev[6],40)))
all_evy=data.frame(all_years,all_evt)
all_evy$yev=paste(all_evy$year,all_evy$event,sep=" ")
event_tag$yev=paste(event_tag$year,event_tag$Events,sep=" ")
df_filled <- left_join(all_evy, event_tag, by = "yev")
for (c in 1:length(df_filled)){
  df_filled[which(is.na(df_filled[,c])),c]=0
}
event_tag$avg10=NA



# mylab= "total affected percentage"
# mylab= "total losses percentage"
# mylab= "total fatalities percentage"
# mylab= "percentage of events"


data <- df_filled  %>%
  group_by(year.x, event) %>%
  summarise(n = sum(losses.len)) %>%
  mutate(percentage = n / sum(n)*100)

data$avg10=NA
for (evn in unique(df_filled$event)){
  event_sub=data[which(data$event==evn),]
  event_sub=event_sub[order(event_sub$year.x),]
  event_sub$avg10=tsEvaNanRunningMean(event_sub$percentage,10)
  # plot(event_sub$avg10)
  data$avg10[which(data$event==evn)]=event_sub$avg10
}
ggplot()+
  geom_bar(data=event_tag, aes(x = year, y = losses.len, fill = Events),stat = "identity",color="transparent", alpha=0.9) +
  scale_y_continuous(
    name = "Number of events",
    breaks=seq(0,200,20),
    minor_breaks = seq(0,200,10),
    sec.axis = sec_axis( transform=~.*.5, name="Proportion of events (%)",
                         breaks=seq(0,100,20))
  )+
  geom_line(data=data,aes(x=year.x, y=avg10*2,color=event),lwd=1)+
  geom_point( data=data,aes(x=year.x, y=avg10*2,fill=event),size = 1.5,
              pch = 21, # Type of point that allows us to have both color (border) and fill.
              color = "white",
              stroke = .8 )+
  #geom_area(aes(year, percentage, fill = Events), color = "black") +
  #geom_point(data=taglist,aes(x=factor(Events),y=losses.mean),color="red",size=3,position=position_dodge(.9))+
  #scale_y_continuous(name="log(losses)")+
  # scale_x_discrete(labels=c("1" = "100-200", "2" = "200-500",
  #                           "3" = "500-1000","4" = "1000-10 000",
  #                           "5" = "10 000-100 000","6" = ">100 000"),name="Catchment Area (km2)")+
  scale_color_manual(name="Event archetypes",
                     values=color_events) +
  scale_fill_manual(name="Event archetypes",
                    values=color_events,drop=FALSE) +
  labs(x="Year",y=mylab) +
  #geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/change_in_time_F2.jpg"), width=40, height=15, units=c("cm"),dpi=300) 



color_events=c("single-flood"="gray","wet-sequence"="dodgerblue4","heat-flood"="darksalmon","drought-flood"="orange","cold-flood"="lightblue","windstorm-flood"="seagreen")


##Plot of seasonality for each event ----

event_class3$Type[which(event_class3$Type=="River/Coastal")]="River"
event_class3$yday=yday(as.Date(event_class3$Start.date))
event_class3$theta=event_class3$yday*(2*pi/365.25)

ydays = bind_rows(event_class3, event_class3, event_class3)
ydays$bday = ydays$yday + rep(c(0,365,365*2),each=length(event_class3$yday))

dtf=event_class3[which(event_class3$event_name=="drought-flood"),]
seasondtf=seasony(dtf$yday)
length(dtf$eventID)

sing=event_class3[which(event_class3$event_name=="single-flood"),]
seasonsing=seasony(sing$yday)
length(sing$eventID)

fsq=event_class3[which(event_class3$event_name=="wet-sequence"),]
seasonfsq=seasony(fsq$yday)
length(fsq$eventID)

htf=event_class3[which(event_class3$event_name=="heat-flood"),]
seasonhtf=seasony(htf$yday)
length(htf$eventID)

ctf=event_class3[which(event_class3$event_name=="cold-flood"),]
seasonctf=seasony(ctf$yday)
length(ctf$eventID)

seasonall=seasony(event_class3$yday)

newday=ydays[c((length(event_class3$yday)+1):(2*length(event_class3$yday))),]


unique(ydays$Type)
ydays$Type[which(ydays$Type=="River/Coastal")]="River"
color_type=c("River/Coastal"="purple","River"="slateblue","Flash"="yellow3")

ydays$event_name <- factor(ydays$event_name, 
                         levels = c("single-flood","wet-sequence","heat-flood","drought-flood","cold-flood","windstorm-flood"),
                         ordered = TRUE)

ptn=direction_labeller(ydays$bday)

###[Plot] - seasonnality of different compound events ----
ggplot(ydays, aes(y=factor(event_name), x=(bday),color=Type)) +
  # geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=event_name),linewidth=0.8,outlier.alpha = 0.4)+
  geom_quasirandom(alpha=0.5, size=3, bandwidth = 0.1)+
  coord_cartesian(xlim=c(366,731))+
  scale_y_discrete(name="MH events")+
  scale_size(range = c(1, 10), trans="sqrt",name= "losses (million euros)",
             breaks=c(100000,1000000,10000000,100000000,1000000000,10000000000), labels=c("0.1","1", "10", "100", "1000","10 000"))+
  # scale_size(range = c(0.08, 0.4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
  #                                                                      sep = " ")),
  #            breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))
  scale_x_continuous(name="day of the year",labels=direction_labeller, breaks=seq(732,30.5,-30.5),expand = c(0, 0))+
  scale_color_manual(
    values=color_type) +
  guides(colour = guide_legend(override.aes = list(size = 10)))+
  #geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=18, face="bold"),
        axis.text = element_text(size=14),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        legend.position = "top",
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.x = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/compound_events_seasonality_informed5.png"), width=30, height=20, units=c("cm"),dpi=300) 



###[Plot] - Piechart or event type per compound ------

event_pie<-aggregate(list(losses=event_class3$Losses..2020.euro.),
                     by = list(Events=event_class3$event_name,Type=event_class3$Type),
                     FUN = function(x) c(len=length(x)))
event_pie <- do.call(data.frame, event_pie)

event_p2=event_pie[which(event_pie$Events=="cold-flood"),]
event_p3=event_pie[which(event_pie$Events=="drought-flood"),]
event_p4=event_pie[which(event_pie$Events=="heat-flood"),]
event_p5=event_pie[which(event_pie$Events=="wet-sequence"),]
event_p1=event_pie[which(event_pie$Events=="windstorm-flood"),]
event_p6=event_pie[which(event_pie$Events=="single-flood"),]

p1 <- ggplot(event_p1, aes(x="", y=losses, fill=Type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual(
    values=color_type) +
  labs(title = "windstorm-flood")+
  theme_void()+
  theme(legend.position = "none")

p2 <- ggplot(event_p2, aes(x="", y=losses, fill=Type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual(
    values=color_type) +
  labs(title = "cold-flood")+
  theme_void()+
  theme(legend.position = "none")

p3 <- ggplot(event_p3, aes(x="", y=losses, fill=Type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual(
    values=color_type) +
  labs(title = "drought-flood")+
  theme_void()+
  theme(legend.position = "none")

p4 <- ggplot(event_p4, aes(x="", y=losses, fill=Type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  labs(title = "heat-flood")+
  scale_fill_manual(
    values=color_type) +
  theme_void()+
  theme(legend.position = "none")

p5 <- ggplot(event_p5, aes(x="", y=losses, fill=Type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  labs(title = "wet-sequence")+
  scale_fill_manual(
    values=color_type) +
  theme_void()+
  theme(legend.position = "none")

p6 <- ggplot(event_p6, aes(x="", y=losses, fill=Type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  labs(title = "single-flood")+
  scale_fill_manual(
    values=color_type) +
  theme_void()+
  theme(legend.position = "none")



# Create the grid of pie charts
p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 1)


ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/pie_evtypes2.jpg"), width=10, height=30, units=c("cm"),dpi=300) 


### Plot map of Hanze event frequency by NUTS3 ----



### Plot map of event type frequency by NUTS3 ----
#Merge with NUTS shapefile
event_MHplot=right_join(event_MHRn,NUTS3,by=c("NUT"="NUTS_ID"))

event_MHplot$bin="0"
event_MHplot$bin[which(!is.na(event_MHplot$mh))]="1-2"
event_MHplot$bin[which(event_MHplot$mh>2)]="2-5"
event_MHplot$bin[which(event_MHplot$mh>5)]="5-10"
event_MHplot$bin[which(event_MHplot$mh>10)]="10-20"
event_MHplot$bin[which(event_MHplot$mh>20)]=">20"

#Plot

labels=c("0","1-2","2-5","5-10","10-20",">20")
legend2="# multi-hazard events"
palet=c(hcl.colors(6, palette = "YlOrRd", alpha = NULL, rev = T, fixup = TRUE))
paletf=c(hcl.colors(11, palette = "Purples", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white",color="darkgrey",size=0.5)+
  geom_sf(data=event_MHplot,aes(fill=bin,geometry=geometry),alpha=0.9,color="transparent")+
  scale_fill_manual(
    values=palet, breaks=rev(labels), name=legend2)   +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_color_gradientn(
  #   colors=palet,
  #   breaks=br,
  #   oob = scales::squish,na.value="transparent", name=legend2)   +
  labs(x="Longitude", y = "Latitude")+
  #guides(fill = guide_coloursteps(barwidth = 1.5, barheight = 14,nbins=5))+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=osize),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/MH_events_8120_HANZE.jpg"), width=30, height=20, units=c("cm"),dpi=300) 



#Merge with NUTS shapefile
event_MHplot=right_join(event_Nt2,NUTS3,by=c("NUT"="NUTS_ID"))

event_MHplot$bin=as.character(event_MHplot$ev)
event_MHplot$bin[which(is.na(event_MHplot$bin))]="0"

###[Plot]
palet=c(hcl.colors(7, palette = "YlOrRd", alpha = NULL, rev = F, fixup = TRUE))
labels=c("0","1","2","3","4","5","6")
legend2="# Event archetypes"
ggplot(basemap) +
  geom_sf(fill="white",color="darkgrey",size=0.5)+
  geom_sf(data=event_MHplot,aes(fill=bin,geometry=geometry),alpha=0.9,color="transparent")+
  scale_fill_manual(
    values=palet, breaks=rev(labels), name=legend2)   +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_color_gradientn(
  #   colors=palet,
  #   breaks=br,
  #   oob = scales::squish,na.value="transparent", name=legend2)   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_legend(title.position = "top",nrow = 1,
                             reverse=T, title.hjust=0.5, 
                             barwidth = 5, barheight = 1.5)) +
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=tsize),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",  # Stack legends vertically
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"),
        legend.key.spacing.x = unit(0, "cm"),
        legend.text.position = "bottom",
        legend.margin = margin(t = 0, unit = "cm") )# Reduce the spacing between keys

  

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/MH_events_types_8120_HANZE2.jpg"), width=20, height=30, units=c("cm"),dpi=400) 




#Aggregate HANZE events by NUTS

HANZE_NA<-aggregate(list(losses= Hanze_flood_BN2$Losses..2020.euro.,fat= Hanze_flood_BN2$Fatalities,af=Hanze_flood_BN2$Persons.affected),
                     by = list(NUTS= Hanze_flood_BN2$NUTl),
                     FUN = function(x) c(mean=mean(x,na.rm=T),len=length(x),sum=sum(x,na.rm=T)))
HANZE_NA <- do.call(data.frame, HANZE_NA)

#Merge with NUTS shapefile
Hanze_NAp=right_join(HANZE_NA,NUTS3,by=c("NUTS"="NUTS_ID"))

Hanze_NAp$bin="0"
Hanze_NAp$bin[which(!is.na(Hanze_NAp$fat.len))]="1-2"
Hanze_NAp$bin[which(Hanze_NAp$fat.len>2)]="2-5"
Hanze_NAp$bin[which(Hanze_NAp$fat.len>5)]="5-10"
Hanze_NAp$bin[which(Hanze_NAp$fat.len>10)]="10-20"
Hanze_NAp$bin[which(Hanze_NAp$fat.len>20)]=">20"

###[Plot] - total number of events per NUTS3 ----
labels=c("0","1-2","2-5","5-10","10-20",">20")
legend2="# events"
palet=c(hcl.colors(6, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE))

paletf=c(hcl.colors(11, palette = "Purples", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white",color="darkgrey",size=0.5)+
  geom_sf(data=Hanze_NAp,aes(fill=bin,geometry=geometry),alpha=0.9,color="transparent")+
  scale_fill_manual(
    values=palet, breaks=rev(labels), name=legend2)   +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_color_gradientn(
  #   colors=palet,
  #   breaks=br,
  #   oob = scales::squish,na.value="transparent", name=legend2)   +
  labs(x="Longitude", y = "Latitude")+
  #guides(fill = guide_coloursteps(barwidth = 1.5, barheight = 14,nbins=5))+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=osize),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/All_floods_8120_HANZE.jpg"), width=30, height=20, units=c("cm"),dpi=300) 



#losses plot:
Hanze_NAp$losses.sumM=Hanze_NAp$losses.sum/1e6

paletf=c(hcl.colors(11, palette = "Purples", alpha = NULL, rev = T, fixup = TRUE))
legend2="Total losses \n(million Euros 2020)"
palet=c(hcl.colors(11, palette = "viridis", alpha = NULL, rev = F, fixup = TRUE))

limi=c(1e0,1e4)
Hanze_NAp$losses.sum[which(Hanze_NAp$losses.sum==0)]=NA

ggplot(basemap) +
  geom_sf(fill="white",color="darkgrey",size=0.5)+
  geom_sf(data=Hanze_NAp,aes(fill=losses.sumM,geometry=geometry),alpha=0.9,color="transparent")+
  scale_fill_gradientn(
    colors=palet,
    breaks=c(1e0,1e1,1e2,1e3,1e4), labels=c(0,10,100,1000,10000), limits=limi,trans="log",
    oob = scales::squish,na.value=palet[1], name=legend2)   +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_color_gradientn(
  #   colors=palet,
  #   breaks=br,
  #   oob = scales::squish,na.value="transparent", name=legend2)   +
  labs(x="Longitude", y = "Latitude")+
  #guides(fill = guide_coloursteps(barwidth = 1.5, barheight = 14,nbins=5))+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=osize),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))


####[Plot] - Dominant event archetype ----

event_MHplot=inner_join(event_Nt2,NUTS3,by=c("NUT"="NUTS_ID"))

color_events=c("single-flood"="gray","wet-sequence"="dodgerblue4","heat-flood"="darksalmon","drought-flood"="orange","cold-flood"="lightblue","windstorm-flood"="seagreen")

legend2="Dominant event type"
# palet=c(hcl.colors(7, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE))

paletf=c(hcl.colors(11, palette = "Purples", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white",color="darkgrey",size=0.5)+
  geom_sf(data=event_MHplot,aes(fill=dominant,geometry=geometry),alpha=0.9,color="transparent")+
  scale_fill_manual(
    values=color_events, name=legend2)   +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  # scale_color_gradientn(
  #   colors=palet,
  #   breaks=br,
  #   oob = scales::squish,na.value="transparent", name=legend2)   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_legend(title.position = "top",nrow = 1,
                             reverse=T, title.hjust=0.5, 
                             barwidth = 6, barheight = 1.5)) +
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=tsize),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",  # Stack legends vertically
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"),
        legend.key.spacing.x = unit(0, "cm"),
        legend.text.position = "bottom",
        legend.margin = margin(t = 0, unit = "cm") )# Reduce the spacing between keys

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/dominant_events_types_8120_HANZE2.jpg"), width=20, height=30, units=c("cm"),dpi=400) 


length(which(event_MHplot$dominant=="single-flood"))/length(event_MHplot$dominant[which(!is.na(event_MHplot$dominant))])


## Impact class analysis ----
event_classD=event_class3
event_classD$damage=event_classD$Losses..2020.euro.
event_classD$damage[which(is.na(event_classD$damage))]=0
event_classD$category=event_classD$event_name

oo=ecdf(event_classD$damage)
ld=length(event_classD$eventID)
dama=data.frame(event=event_classD$eventID,damage=event_classD$damage)
dama=dama[order(-dama$damage),]
dama$id=seq(1, ld)
plot(dama$damage)
total_damage=sum(dama$damage)
cumulative_percentages <- (cumsum(dama$damage) / total_damage) * 100
cumulative_length <- seq(1, ld)
plot(cumulative_length,cumulative_percentages)
plot(dama$damage,cumulative_percentages)
dama$cumul=cumulative_percentages

damage_data <- data.frame(
  Event_ID = event_classD$eventID,
  Damage_Amount = event_classD$damage,
  Category = event_classD$event_name# Replace with your actual categories
)

Dam=event_classD$damage
# Print summary of the fit


total_damage <- sum(damage_data$Damage_Amount)

cumulative_percentages <- numeric(nrow(damage_data))
proportions <- matrix(nrow = nrow(damage_data), ncol = length(unique(damage_data$Category)))

cumulative_data <- data.frame(
  Cumulative_Percentage = numeric(),
  Category = numeric(),
  Proportion = numeric(),
  Proportion_ev = numeric(),
  Number_of_Events = numeric()
)

damage_data=damage_data[order(-damage_data$Damage_Amount),]

###Plots on cumulative frequency ----

for (i in 1:nrow(damage_data)) {
  # Calculate the cumulative damage of the top i events
  cumulative_damage <- sum(damage_data$Damage_Amount[1:i])
  
  # Calculate the cumulative percentage
  cumulative_percentage <- (cumulative_damage / total_damage) * 100
  # Calculate the proportion of each category in the cumulative damage
  cumulative_categories <- damage_data[1:i, ]
  for (j in unique(damage_data$Category)) {
    category_damage <- sum(cumulative_categories$Damage_Amount[cumulative_categories$Category == j])
    proportion <- category_damage / cumulative_damage
    number_of_events <- sum(cumulative_categories$Category == j)
    proportion2 <- number_of_events / length(cumulative_categories$Category)
    
    # Append the data to the cumulative_data data frame
    cumulative_data <- rbind(cumulative_data, data.frame(
      Cumulative_Percentage = cumulative_percentage,
      length=i,
      Category = j,
      Proportion = proportion,
      Proportion_ev = proportion2,
      Number_of_Events = number_of_events
    ))
  }
  
}

# Print the cumulative data
print(cumulative_data)
ld=length(damage_data$Event_ID)
#cumulative_data$Cumulative_Percentage=100-cumulative_data$Cumulative_Percentage
#use approx to know to which value of cumulative percentage number of events is correponding to
localoca=approx(cumulative_data$length,cumulative_data$Cumulative_Percentage,xout=c(500,200,100,50,10,1))
ggplot(cumulative_data, aes(x = Cumulative_Percentage, y = Proportion_ev, color = factor(Category))) +
  geom_line(size = 1) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(
    name = "Cumulative Percentage of Total Damage",
    breaks = c(0,20,40,60,80,100),
    sec.axis = sec_axis(~.*1,name = "Cumulative Number of Events", 
                        breaks = localoca$y,
                        labels = localoca$x)) +
  labs(y = "Proportion of Events",
       color = "Category") +
  scale_color_manual(
    values=color_events) +
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.title.x.top = element_text(hjust = 0.5),
        axis.title.x = element_text(hjust = 0.5),
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))




for (i in 1:nrow(damage_data)) {
  cumulative_damage <- sum(damage_data %>% arrange(desc(Damage_Amount)) %>% slice(1:i) %>% pull(Damage_Amount))
  cumulative_percentage <- (cumulative_damage / total_damage) * 100
  cumulative_percentages[i] <- cumulative_percentage
  
  cumulative_categories <- damage_data %>% arrange(desc(Damage_Amount)) %>% slice(1:i)
  proportions[i, ] <- sapply(unique(damage_data$Category), function(j) {
    sum(cumulative_categories$Damage_Amount[cumulative_categories$Category == j]) / cumulative_damage
  })
}
length(cumulative_length)
# Cairo::Cairo(
#   20, #length
#   15, #width
#   file = "D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/Hanze_Cumlosses.png",
#   type = "png", #tiff
#   bg = "transparent", #white or transparent depending on your requirement 
#   dpi = 200,
#   units = "cm" #you can change to pixels etc 
# )

plot(cumulative_length,cumulative_percentages, ylab="cumulative loss percentages",lwd=2.5,
     xlab="cumulative numver of events", type="n", ylim=c(0,100),main="loss distribution in Hanze events (n=689)")
grid(nx = NULL, ny = NULL,
     lty = 2, col = "gray", lwd = 1)
lines(cumulative_length,cumulative_percentages, lwd=2.5)
dev.off()
print(data.frame(Cumulative_Percentage = cumulative_percentages))
colnames(proportions)=unique(damage_data$Category)
damage_data=data.frame(cumulative_length,cumulative_percentages,proportions)

cziz=data.frame(cumulative_length,cumulative_percentages)


plot(100-damage_data$cumulative_percentage,damage_data$cold.flood)
# Calculate percentiles
event_classD2=event_classD[which(event_classD$damage>0),]
percentiles <- c(0,seq(1, 99, by = 1))
damage_percentiles <- quantile(event_classD2$damage, probs = percentiles / 100)
damage_percentiles[1]<-0



# Calculate proportions


#attribute bins to event_classD
event_classD$bins=0
event_classD$bins[which(event_classD$Losses..2020.euro.>=damage_percentiles[1] & event_classD$Losses..2020.euro.<damage_percentiles[51])]=1
event_classD$bins[which(event_classD$Losses..2020.euro.>=damage_percentiles[51] & event_classD$Losses..2020.euro.<damage_percentiles[91])]=2
event_classD$bins[which(event_classD$Losses..2020.euro.>=damage_percentiles[91] & event_classD$Losses..2020.euro.<damage_percentiles[100])]=3
event_classD$bins[which(event_classD$Losses..2020.euro.>=damage_percentiles[100])]=4

proportions_data <- calculate_proportions(df=event_classD, percentiles, damage_percentiles)

# Plot
ggplot(proportions_data, aes(x = percentile/100, y = proportion, color = category, group = category)) +
  geom_line(size = 1) +
  geom_point(alpha = 0.5) +
  labs(title = "Proportion of Events Above Percentiles by Category Among All Events",
       x = "Loss Percentile",
       y = "Proportion Above Percentile") +
  scale_color_manual(
    values=color_events) +
  coord_trans(x = "exp") +
  #geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/compound_events_LossP",wd,".jpg"), width=30, height=20, units=c("cm"),dpi=300) 

#bins of top 1%, 2/10%, 10/50%, remaining 50%

proportions_data$bins=4
proportions_data$bins[which(proportions_data$percentile<99 & proportions_data$percentile>=90)]=3
proportions_data$bins[which(proportions_data$percentile<90 & proportions_data$percentile>=50)]=2
proportions_data$bins[which( proportions_data$percentile<50)]=1
proportions_data$bins[which( proportions_data$percentile<1)]=0
proportions_data[which(proportions_data$percentile<99 & proportions_data$percentile>=90),]

#plot by beans

frequencies <- proportions_data %>%
  group_by(bins, category) %>%
  summarise(mean = mean(proportion), .groups = 'drop')

frequenciesCor <- event_classD %>%
  group_by(bins,event_name) %>%
  summarise(mx =length(nh), .groups = 'drop')

frequenciesCB <- event_classD %>%
  group_by(bins) %>%
  summarise(mx =length(eventID), .groups = 'drop')

sum(frequenciesCB$mx)
mb=match(frequenciesCor$bins,frequenciesCB$bins)
frequenciesCor$lbin=frequenciesCB$mx[mb]

frequenciesCor$prop=frequenciesCor$mx/frequenciesCor$lbin
frequencies2 <- proportions_data %>%
  group_by(bins) %>%
  summarise(mx = max(damage_perc),mx2=max(length), .groups = 'drop')

frequencies2$mx3=frequencies2$mx2

for (ic in 3:1){
  frequencies2$mx3[ic]=frequencies2$mx3[ic]-frequencies2$mx2[ic+1]
}
impact_leg=c("no recorded impact","lower 50%","50% - 10%","10% - 1%"," top 1%")
# Plotting
library(ggplot2)

ggplot() +
  geom_bar(data=frequenciesCor, aes(x = bins, y = prop, fill = event_name),stat = "identity",color="black") +
  geom_line(data=frequencies2,aes(x=bins,y=mx),lwd=1)+
  geom_point(data=frequencies2,aes(x=bins,y=mx),pch=21,fill="white",stroke=1,size=2)+
  labs(
    x = "Impact class",
    y = "Proportion",
    fill = "Event types") +
  scale_x_continuous(breaks=c(0,1,2,3,4),labels=impact_leg)+
  scale_y_continuous(
    name = "Proportion of events",
    breaks=seq(0,1,.2),
    minor_breaks = seq(0,1,.1),
    sec.axis = sec_axis( transform=~.*100, name="Cumulatice percentage of total losses",
                         breaks=seq(0,100,20))
  )+
  scale_fill_manual(
    values=color_events) +
  geom_text(data=frequenciesCor, aes(x=bins, y=1.05, label=lbin), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=14),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        #panel.grid.major = element_line(colour = "grey60"),
        #panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/compound_events_Impact_classesF2.jpg"), width=30, height=20, units=c("cm"),dpi=300) 




#RDH VI map

RDH_VI= read.csv("D:/tilloal/Documents/01_Projects/RiskDynamics/Data/multihz_data/rdh_vi.csv")
RDH_2020=RDH_VI[which(RDH_VI$Time=="2020"),]

RDH_plot=inner_join(RDH_2020,NUTS3,by=c("NID"="NUTS_ID"))

legend2="RDH VI (2020)"
palet=c(hcl.colors(8, palette = "viridis", alpha = NULL, rev = F, fixup = TRUE))
paletf=c(hcl.colors(11, palette = "Purples", alpha = NULL, rev = T, fixup = TRUE))

ggplot(basemap) +
  geom_sf(fill="white",color="darkgrey",size=0.5)+
  geom_sf(data=RDH_plot,aes(fill=VI,geometry=geometry),alpha=0.9,color="transparent")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet,
    breaks=c(1,2,3,4,5,6),
    oob = scales::squish,na.value="white", name=legend2)   +
  # scale_color_gradientn(
  #   colors=palet,
  #   breaks=br,
  #   oob = scales::squish,na.value="transparent", name=legend2)   +
  labs(x="Longitude", y = "Latitude")+
  #guides(fill = guide_coloursteps(barwidth = 1.5, barheight = 14,nbins=5))+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=osize),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))

ggsave(paste0("D:/tilloal/Documents/01_Projects/RiskDynamics/Plots/RDH_VI_2020.jpg"), width=30, height=20, units=c("cm"),dpi=300) 


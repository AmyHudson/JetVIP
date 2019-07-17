#JetHypothesisFig_JA_scpdsi.R

# Jet Stream Hypothesis Figure: I envision April May Jet stream as faint lines over the northern hemisphere background, with extreme poleward (10) jets bolded and extreme equatorward bolded. 

#Load Jet stream


# Read in JA Jet stream 
jet <- read.table("JetIndices.txt", header = TRUE)
jet <- jet[jet$YEAR>=1981,]
indices_names <- colnames(jet)
which(indices_names == "JA_Reg1") #jet[,19]

library(maps)
library(mapdata)
library(maptools)  #for shapefiles
library(scales)  #for transparency
library(ggplot2)
library(devtools)
library(dplyr)
library(stringr)
library(ggmap) #need to cite ggmap
library(RColorBrewer)
library(raster)


latmax <- 70
latmin <- 24
lonmax <- 180
lonmin <- -180

##################################################################################

# Now composite relationship with temperature field for subsetted region.

setwd("/Volumes/AOP-NEON1.4")

library(ncdf4)
library(fields)
library(Hmisc)
library(mapdata)

#read in scPDSI field
data_time <- 1901:2017
year <- matrix(rep(1901:2017,12),nrow = length(1901:2017),ncol = 12)
year <- as.vector(t(year))
month <- rep(1:12,117)#repeat 1:12 116 times
time <- cbind(year,month)

a <- nc_open("scPDSI.cru_ts3.26early.bams2018.GLOBAL.1901.2017.nc")

xdim <- round(a$dim[[1]]$vals, digits = 5)# lon
ydim <- round(a$dim[[2]]$vals, digits = 5) # lat 
zdim <- a$dim[[3]]$vals # time

xs <- which(xdim == -179.75)
ys <-  which(ydim == 69.25) #89.8:-89.8 
# this is opposite of the temperature data set, which only comes up because we need to flip the temperature raster: rotate3 <- raster(rho1[nrow(rho1):1,]) here, we just assign rotate3 <- raster(rho1)
#zs <- which(zdim == 34348)

pdsi <- a$var[[1]]

ptm <- proc.time()   
#Clim_var.nc<- ncvar_get(a, pdsi,start = c(xs,ys,zs), count = c(114, 92,276))  
Clim_var.nc<- ncvar_get(a, pdsi,start = c(xs,ys,1), count = c((which(xdim == 179.75)-which(xdim == -179.75)+1), (which(ydim == 24.25)-which(ydim == 69.25)+1),1404))    
proc.time() - ptm
dim(Clim_var.nc)[1] #720
dim(Clim_var.nc)[2] #91
dim(Clim_var.nc)[3] #1404

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]
latitude <- ydim[c(ys:(ys+(which(ydim == 24.25)-which(ydim == 69.25))))]

mon1 <- Clim_var.nc[,,month == 7] 
mon2 <- Clim_var.nc[,,month == 8]
mon <- (mon1 + mon2) / 2

pdsi1 <- aperm(mon,c(3,2,1)) #reorder with time variable in front, lat, lon
#time is now in years

###########################################################################################
# Grab 10 north/south jet streams JA Region 1: 20-54E
indices_names[19]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

all <- as.data.frame(cbind(jet$YEAR,jet[,19]))
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 54 & longitude >=20)]
#crop to 1981:2012
#crop to 24 to 74 AM_1
longitude <- longitude[which(longitude <= 54 & longitude >=20)]
#latitude <- latitude[which(latitude <= 69.25 & latitude >=24.25)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) #reshape 2D w time, lat*lon
pdsi2df <- as.data.frame(pdsi2)

year <- 1981:2012

#m1 <- pdsi2df %>% filter(year %in% as.character(allN))
m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

my.t.test.p.value <- function(...) { 
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "greater",na.omit = T) #greater for jet poleward
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j1 <- rotate3

####################################################################################################
# JA2
indices_names[20]

all <- as.data.frame(cbind(jet$YEAR,jet[,20]))
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 86 & longitude >=54)]
longitude <- longitude[which(longitude <= 86 & longitude >=54)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

#m1 <- pdsi2df %>% filter(year %in% as.character(allN))
m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "greater",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j2 <- rotate3
map("world", xlim=c(54,86),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j2,add = T)
###################################################################################################
# JA3
indices_names[21]
all <- as.data.frame(cbind(jet$YEAR,jet[,21])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 180 & longitude >=86)]
longitude <- longitude[which(longitude <= 180 & longitude >=86)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

#m1 <- pdsi2df %>% filter(year %in% as.character(allN))
m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "greater",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j3 <- rotate3
#map("world", xlim=c(86,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

#plot(j3,add = T)
###################################################################################################
# JA3 
q <- 21
lonmin1 <- -180
lonmax1 <- -160

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

#m1 <- pdsi2df %>% filter(year %in% as.character(allN))
m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "greater",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j3b <- rotate3
############################################################################
###################################################################################################
# JA4
q <- 22
lonmin1 <- -160
lonmax1 <- -104

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

#m1 <- pdsi2df %>% filter(year %in% as.character(allN))
m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "greater",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j4 <- rotate3

###################################################################################################
# JA5
q <- 23
lonmin1 <- -104
lonmax1 <- -58

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

#m1 <- pdsi2df %>% filter(year %in% as.character(allN))
m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "greater",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

j5 <- rotate3

###################################################################################################
# JA6 
q <- 24
lonmin1 <- -58
lonmax1 <- -16

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

#m1 <- pdsi2df %>% filter(year %in% as.character(allN))
m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "greater",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

j6 <- rotate3

###################################################################################################
# JA7
q <- 25
lonmin1 <- -16
lonmax1 <- 8

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

#m1 <- pdsi2df %>% filter(year %in% as.character(allN))
m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "greater",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j7 <- rotate3

###################################################################################################

##########################################################################################################################
t1 <- merge(j1,j2,tolerance = 0.5)
#map("world", xlim=c(20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
plot(j3,add = T,legend = F)
plot(t1,add = T,legend = F)


#nothing in j6
t <- merge(j4,j5,j6, tolerance = 0.5) 
map("world", xlim=c(-180,8),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
#t2 <- merge(j4,j5, tolerance = 0.5) 
plot(t,add = T, legend = F)
plot(j7,add = T, legend = F)

# Plot all 

par(mai=c(0,0,0,0))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 


for (i in 1:length(jet$JA_Reg1)){
  lines(20:54,rep(jet$JA_Reg1[i],length(20:54)), col = "grey")
}
for (i in 1:length(jet$JA_Reg2)){
  lines(54:86,rep(jet$JA_Reg2[i],length(54:86)), col = "grey")
}
for (i in 1:length(jet$JA_Reg3)){
  lines(86:180,rep(jet$JA_Reg3[i],length(86:180)), col = "grey")
}
for (i in 1:length(jet$JA_Reg3)){
  lines(-180:-160,rep(jet$JA_Reg3[i],length(-180:-160)), col = "grey")
}
for (i in 1:length(jet$JA_Reg4)){
  lines(-160:-104,rep(jet$JA_Reg4[i],length(-160:-104)), col = "grey")
}
for (i in 1:length(jet$JA_Reg5)){
  lines(-104:-58,rep(jet$JA_Reg5[i],length(-104:-58)), col = "grey")
}
for (i in 1:length(jet$JA_Reg6)){
  lines(-58:-16,rep(jet$JA_Reg6[i],length(-58:-16)), col = "grey")
}
for (i in 1:length(jet$JA_Reg7)){
  lines(-16:8,rep(jet$JA_Reg7[i],length(-16:8)), col = "grey")
}

map.axes() 


rb <- rev(brewer.pal(n = 7, name = 'RdBu'))
rb <- c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7","#F7F7F7", "#FDDBC7", "#EF8A62","#B2182B")
bg <- brewer.pal(n = 7, name = 'BrBG')

par(mai=c(0,0,0,0))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg1)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(20:54,rep(S[i],length(20:54)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg2)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(54:86,rep(S[i],length(54:86)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg3)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(86:180,rep(S[i],length(86:180)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg3)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-180:-160,rep(S[i],length(-180:-160)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg4)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-160:-104,rep(S[i],length(-160:-104)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg5)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-104:-58,rep(S[i],length(-104:-58)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg6)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-58:-16,rep(S[i],length(-58:-16)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg7)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-16:8,rep(S[i],length(-16:8)), col = "grey")
}


abline(v = 20)
abline(v = 54)
abline(v = 86)
abline(v = -160)
abline(v = -104)
abline(v = -58)
abline(v = -16)
abline(v = 8)

map.axes()

plot(t,main="",col =bg,breaks = c(-2.5,-1.5,-1,-0.5,0.5,1,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j7,main="",col =bg,breaks = c(-2.5,-1.5,-1,-0.5,0.5,1,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)

#plot(t1,main="",col =bg,breaks = c(-2,-1.5,-1,-0.5,0.5,1,1.5,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)
#plot(j3,main="",col =bg,breaks = c(-2,-1.5,-1,-0.5,0.5,1,1.5,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)

plot(t,main="",col =bg,breaks = c(-2,-1.5,-1,-0.5,0.5,1,1.5,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j7,main="",col =bg,breaks = c(-2,-1.5,-1,-0.5,0.5,1,1.5,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)

##############################################################################################################################

## NOW POLEWARD!
# Grab 10 north/south jet streams JA Region 1: 20-54E
indices_names[19]
longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

all <- as.data.frame(cbind(jet$YEAR,jet[,19]))
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 54 & longitude >=20)]
#crop to 1981:2012
#crop to 24 to 74 AM_1
longitude <- longitude[which(longitude <= 54 & longitude >=20)]
#latitude <- latitude[which(latitude <= 69.25 & latitude >=24.25)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) #reshape 2D w time, lat*lon
pdsi2df <- as.data.frame(pdsi2)

year <- 1981:2012

m1 <- pdsi2df %>% filter(year %in% as.character(allN))
#m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

my.t.test.p.value <- function(...) { 
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "less",na.omit = T) #less for jet poleward
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j1 <- rotate3

####################################################################################################
# JA2
indices_names[20]

all <- as.data.frame(cbind(jet$YEAR,jet[,20]))
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 86 & longitude >=54)]
longitude <- longitude[which(longitude <= 86 & longitude >=54)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

m1 <- pdsi2df %>% filter(year %in% as.character(allN))
#m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "less",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j2 <- rotate3
map("world", xlim=c(54,86),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j2,add = T)
###################################################################################################
# JA3
indices_names[21]
all <- as.data.frame(cbind(jet$YEAR,jet[,21])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 180 & longitude >=86)]
longitude <- longitude[which(longitude <= 180 & longitude >=86)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

m1 <- pdsi2df %>% filter(year %in% as.character(allN))
#m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "less",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j3 <- rotate3
#map("world", xlim=c(86,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

#plot(j3,add = T)
###################################################################################################
# JA3 
q <- 21
lonmin1 <- -180
lonmax1 <- -160

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

m1 <- pdsi2df %>% filter(year %in% as.character(allN))
#m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "less",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j3b <- rotate3
############################################################################
###################################################################################################
# JA4
q <- 22
lonmin1 <- -160
lonmax1 <- -104

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

m1 <- pdsi2df %>% filter(year %in% as.character(allN))
#m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "less",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j4 <- rotate3

###################################################################################################
# JA5
q <- 23
lonmin1 <- -104
lonmax1 <- -58

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

m1 <- pdsi2df %>% filter(year %in% as.character(allN))
#m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "less",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

j5 <- rotate3

###################################################################################################
# JA6 
q <- 24
lonmin1 <- -58
lonmax1 <- -16

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

m1 <- pdsi2df %>% filter(year %in% as.character(allN))
#m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "less",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude#-360 #convert latitude of 210 to 300 to -150 to  -60
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

j6 <- rotate3

###################################################################################################
# JA7
q <- 25
lonmin1 <- -16
lonmax1 <- 8

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
pdsi2df <- as.data.frame(pdsi2)

m1 <- pdsi2df %>% filter(year %in% as.character(allN))
#m1 <- pdsi2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(pdsi2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],pdsi2df[,i], alternative = "less",na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- diff[i]
  }
}

rho1 <- as.matrix(rho1)
dim(rho1) <- c(length(latitude),length(longitude))
rho1 <- as.data.frame(rho1)
colnames(rho1) <- longitude
rownames(rho1) <- latitude

rho1 <- as.matrix(rho1)
rotate3 <- raster(rho1) 
#rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j7 <- rotate3

###################################################################################################

##########################################################################################################################
t1 <- merge(j7,j1,j2,tolerance = 0.5)
#map("world", xlim=c(20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(-20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
plot(t1,add = T,legend = F)

t2 <- merge(j3,j3b,tolerance = 0.5)
plot(j3,add = T,legend = F)

#nothing in j6
t <- merge(j4,j5,j6, tolerance = 0.5) 
map("world", xlim=c(-180,8),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
#t2 <- merge(j4,j5, tolerance = 0.5) 
plot(t,add = T, legend = F)
# Plot all 
plot(j3b,add = T, legend = F)

rb <- rev(brewer.pal(n = 7, name = 'RdBu'))
rb <- c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7","#F7F7F7", "#FDDBC7", "#EF8A62","#B2182B")


par(mai=c(0,0,0,1))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg1)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(20:54,rep(S[i],length(20:54)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg2)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(54:86,rep(S[i],length(54:86)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg3)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(86:180,rep(S[i],length(86:180)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg3)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-180:-160,rep(S[i],length(-180:-160)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg4)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-160:-104,rep(S[i],length(-160:-104)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg5)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-104:-58,rep(S[i],length(-104:-58)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg6)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-58:-16,rep(S[i],length(-58:-16)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JA_Reg7)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-16:8,rep(S[i],length(-16:8)), col = "grey")
}


abline(v = 20)
abline(v = 54)
abline(v = 86)
abline(v = -160)
abline(v = -104)
abline(v = -58)
abline(v = -16)
abline(v = 8)

map.axes()

plot(t1,main="",col =bg,breaks = c(-2.5,-1.5,-1,-0.5,0.5,1,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j3,main="",col =bg,breaks = c(-2.5,-1.5,-1,-0.5,0.5,1,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)


plot(t,main="",col =bg,breaks = c(-2.5,-1.5,-1,-0.5,0.5,1,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)


#plot(t1,main="",col =bg,breaks = c(-2,-1.5,-1,-0.5,0.5,1,1.5,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)
#plot(j3,main="",col =bg,breaks = c(-2,-1.5,-1,-0.5,0.5,1,1.5,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)

#plot(t,main="",col =bg,breaks = c(-2,-1.5,-1,-0.5,0.5,1,1.5,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)
#plot(j7,main="",col =bg,breaks = c(-2,-1.5,-1,-0.5,0.5,1,1.5,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)

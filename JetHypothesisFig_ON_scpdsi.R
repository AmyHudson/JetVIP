#JetHypothesisFig_ON_scpdsi.R

# Jet Stream Hypothesis Figure: I envision April May Jet stream as faint lines over the northern hemisphere background, with extreme poleward (10) jets bolded and extreme equatorward bolded. 

#Load Jet stream

setwd("~/Documents/JetVIP")

# Read in ON Jet stream 
jet <- read.table("JetIndices.txt", header = TRUE)
jet <- jet[jet$YEAR>=1981,]
indices_names <- colnames(jet)
which(indices_names == "ON_Reg1") #jet[,19]

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

mon1 <- Clim_var.nc[,,month == 10] 
mon2 <- Clim_var.nc[,,month == 11]
mon <- (mon1 + mon2) / 2

pdsi1 <- aperm(mon,c(3,2,1)) #reorder with time variable in front, lat, lon
#time is now in years

###########################################################################################
# Grab 10 north/south jet streams ON Region 1
q <- 26
lonmin1 <- 10
lonmax1 <- 26

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

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
# ON2
q <- 27
lonmin1 <- 26
lonmax1 <- 56

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

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
map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j2,add = T)
###################################################################################################
# ON3
q <- 28
lonmin1 <- 56
lonmax1 <- 104

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

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
map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j3,add = T)
###################################################################################################
# ON4 
q <- 29
lonmin1 <- 104
lonmax1 <- 148

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

map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j4,add = T)
############################################################################
###################################################################################################
# ON5
q <- 30
lonmin1 <- 148
lonmax1 <- 180

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

j5 <- rotate3

map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j5,add = T)
###################################################################################################
# ON5
q <- 30
lonmin1 <- -180
lonmax1 <- -148

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

j5b <- rotate3

map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j5b,add = T)
###################################################################################################
# ON6  
q <- 31
lonmin1 <- -146
lonmax1 <- -98

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

map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j6,add = T)
###################################################################################################
# ON7
q <- 32
lonmin1 <- -98
lonmax1 <- -66

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

map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j7,add = T)
###################################################################################################
# ON8
q <- 33
lonmin1 <- -66
lonmax1 <- -34

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

j8 <- rotate3

map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j8,add = T)
##########################################################################################################################
# ON9
q <- 34
lonmin1 <- -34
lonmax1 <- -10

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

j9 <- rotate3

map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j9,add = T)
##########################################################################################################################
t1 <- merge(j1,j2,tolerance = 0.5)

t3 <- merge(j3,j4,j5,tolerance = 0.5)
#map("world", xlim=c(0,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(0,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
plot(t3,add = T,legend = F)
plot(j1,add = T,legend = F)
plot(j2,add = T,legend = F)


#nothing in j6
t <- merge(j6,j7, tolerance = 0.5) 
map("world", xlim=c(-180,8),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
t2 <- merge(j8,j9, tolerance = 0.5) 

plot(t,add = T, legend = T)
plot(t2,add = T, legend = T)


# Plot all 

bg <- brewer.pal(n = 5, name = 'BrBG')

par(mai=c(0,0,0,0))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg1)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(10:26,rep(S[i],length(10:26)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg2)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(26:56,rep(S[i],length(26:56)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg3)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(56:104,rep(S[i],length(56:104)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg4)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(104:148,rep(S[i],length(104:148)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg5)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(148:180,rep(S[i],length(148:180)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg5)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-180:-146,rep(S[i],length(-180:-146)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg6)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-146:-98,rep(S[i],length(-146:-98)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg7)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-98:-66,rep(S[i],length(-98:-66)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg8)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-66:-34,rep(S[i],length(-66:-34)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg9)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-34:-10,rep(S[i],length(-34:-10)), col = "grey")
}


abline(v = 10)
abline(v = 26)
abline(v = 56)
abline(v = 104)
abline(v = 148)
abline(v = -146)
abline(v = -98)
abline(v = -66)
abline(v = -34)
abline(v = -10)

map.axes()

plot(j1,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j2,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t3,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)

plot(t,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t2,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)


########################################################################################################################

## NOW POLEWARD!
#ON1
q <- 26
lonmin1 <- 10
lonmax1 <- 26

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

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
# ON2
q <- 27
lonmin1 <- 26
lonmax1 <- 56

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

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
map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j2,add = T)
###################################################################################################
# ON3
q <- 28
lonmin1 <- 56
lonmax1 <- 104

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

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
# ON4 
q <- 29
lonmin1 <- 104
lonmax1 <- 148

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
############################################################################
###################################################################################################
# ON5
q <- 30
lonmin1 <- 148
lonmax1 <- 180

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

j5 <- rotate3

###################################################################################################
# ON5
q <- 30
lonmin1 <- -180
lonmax1 <- -146

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

j5b <- rotate3

###################################################################################################
# ON6 
q <- 31
lonmin1 <- -146
lonmax1 <- -98

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
# ON7
q <- 32
lonmin1 <- -98
lonmax1 <- -66


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
# ON8
q <- 33
lonmin1 <- -66
lonmax1 <- -34


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

j8 <- rotate3
###################################################################################################
# ON9
q <- 34
lonmin1 <- -34
lonmax1 <- -10

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

j9 <- rotate3
##########################################################################################################################

t1 <- merge(j1,j2,tolerance = 0.5)
#map("world", xlim=c(20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(-20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
plot(t1,add = T,legend = F)

t2 <- merge(j3,j4,tolerance = 0.5)
plot(t2,add = T,legend = F)
plot(j5,add = T,legend = F)

#nothing in j6
t <- merge(j6,j7, tolerance = 0.5) 
map("world", xlim=c(-180,8),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
t3 <- merge(j8,j9, tolerance = 0.5) 
plot(t,add = T, legend = F)
# Plot all 
plot(t3,add = T, legend = F)

par(mai=c(0,0,0,1))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg1)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(10:26,rep(S[i],length(10:26)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg2)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(26:56,rep(S[i],length(26:56)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg3)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(56:104,rep(S[i],length(56:104)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg4)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(104:148,rep(S[i],length(104:148)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg5)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(148:180,rep(S[i],length(148:180)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg5)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-180:-146,rep(S[i],length(-180:-146)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg6)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-146:-98,rep(S[i],length(-146:-98)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg7)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-98:-66,rep(S[i],length(-98:-66)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg8)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-66:-34,rep(S[i],length(-66:-34)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg9)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-34:-10,rep(S[i],length(-34:-10)), col = "grey")
}
abline(v = 10)
abline(v = 26)
abline(v = 56)
abline(v = 104)
abline(v = 148)
abline(v = -146)
abline(v = -98)
abline(v = -66)
abline(v = -34)
abline(v = -10)

map.axes()

plot(t1,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t2,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j5,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)

plot(t,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t3,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)

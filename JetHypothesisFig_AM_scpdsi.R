#JetHypothesisFig_AM_scpdsi.R
#scPDSI : a negative value means drier, a positive value means wetter

jet <- read.table("JetIndices.txt", header = TRUE)
jet <- jet[jet$YEAR>=1981,]
indices_names <- colnames(jet)

library(maps)
library(mapdata)
library(maptools)  #for shapefiles
library(scales)  #for transparency
library(ggplot2)
library(devtools)
library(dplyr)
library(stringr)
library(ggmap)
library(RColorBrewer)

library(ncdf4)
library(fields)
library(Hmisc)
library(mapdata)

#read in scPDSI field
data_time <- 1901:2017
#1901:2017 monthly
year <- matrix(rep(1901:2017,12),nrow = length(1901:2017),ncol = 12)
year <- as.vector(t(year))
month <- rep(1:12,117)#repeat 1:12 116 times
#repeat 1901:2017 12 times and make an array for the year vector
#crop to 1994 (jan or feb)
time <- cbind(year,month)
#time <- subset(time,year >=1994)

#don't know how to subset the time variable to 1994
#1117

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

mon1 <- Clim_var.nc[,,month == 4] 
mon2 <- Clim_var.nc[,,month == 5]
mon <- (mon1 + mon2) / 2

pdsi1 <- aperm(mon,c(3,2,1)) #reorder with time variable in front, lat, lon
#time is now in years
####################################################################################################
#AM_1 

# Grab 10 north/south jet streams Region 1: 24 to 74
q <- 10
lonmin1 <- 24
lonmax1 <- 74

all <- as.data.frame(cbind(jet$YEAR,jet[,q]))
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) #reshape 2D w time, lat*lon
pdsi2df <- as.data.frame(pdsi2)

library(dplyr)
#library(tidyverse)

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
#rotate3 <- raster(rho1)
#need to flip rotate 3 on yaxis
#extent(rotate3) <- extent(c(min(longitude-360),max(longitude-360),min(latitude),max(latitude)))
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am1 <- rotate3

##############################################################
###################################################################################################
# NOW AM2

all <- as.data.frame(cbind(jet$YEAR,jet[,11])) #indices_names[11]
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= 116 & longitude >=74)]
longitude <- longitude[which(longitude <= 116 & longitude >=74)]

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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am2 <- rotate3
###################################################################################################
# AM3
indices_names[12]
all <- as.data.frame(cbind(jet$YEAR,jet[,12])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= 152 & longitude >=116)]
longitude <- longitude[which(longitude <= 152 & longitude >=116)]

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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am3 <- rotate3
###################################################################################################
# AM4 #part 1 
q <- 13
lonmin1 <- 152
lonmax1 <- 180

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am4 <- rotate3
############################################################################
###################################################################################################
# AM4 # part 2 
q <- 13
lonmin1 <- -180
lonmax1 <- -150

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am4b <- rotate3

###################################################################################################
# AM5
q <- 14
lonmin1 <- -150
lonmax1 <- -120

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

am5 <- rotate3

###################################################################################################

# AM6 
q <- 15
lonmin1 <- -120
lonmax1 <- -94

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
#rotate3 <- raster(rho1) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

am6 <- rotate3
plot(am6)
###################################################################################################
# AM7
q <- 16
lonmin1 <- -94
lonmax1 <- -56

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am7 <- rotate3

###################################################################################################
# AM8 
q <- 17
lonmin1 <- -10
lonmax1 <- 8

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))


am8 <- rotate3
##########################################################################################################################
t1 <- merge(am1,am2,am3,tolerance = 0.5)
t <- merge(am4,am4b,am5,am6,am7,am8, tolerance = 0.5) 

bg <- brewer.pal(n = 7, name = 'BrBG')

par(mai=c(0,0,0,0))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg1)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(24:74,rep(S[i],length(24:74)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg2)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(74:116,rep(S[i],length(116:74)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg3)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(116:152,rep(S[i],length(116:152)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg4)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(152:180,rep(S[i],length(152:180)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg4)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-180:-150,rep(S[i],length(-180:-150)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg5)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-150:-120,rep(S[i],length(-150:-120)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg6)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-120:-94,rep(S[i],length(-120:-94)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg7)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-94:-56,rep(S[i],length(-94:-56)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg8)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-10:8,rep(S[i],length(-10:8)), col = "grey")
}


abline(v = 24)
abline(v = 74)
abline(v = 116)
abline(v = 152)
abline(v = -150)
abline(v = -120)
abline(v = -94)
abline(v = -56)
abline(v = -10) 
abline(v = 8) 

map.axes()

plot(t1,main="",col =bg,breaks = c(-2,-1.5,-1,-0.5,0.5,1,1.5,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t,main="",col =bg,breaks = c(-2,-1.5,-1,-0.5,0.5,1,1.5,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)

###############################################################################################################
#Now we do the same for the North years
####################################################################################################

#AM_1 

# Grab 10 north/south jet streams Region 1: 24 to 74
q <- 10
lonmin1 <- 24
lonmax1 <- 74

all <- as.data.frame(cbind(jet$YEAR,jet[,q]))
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(pdsi2)<- c(length(1981:2012),length(latitude)*length(longitude)) #reshape 2D w time, lat*lon
pdsi2df <- as.data.frame(pdsi2)

library(dplyr)
#library(tidyverse)

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
#rotate3 <- raster(rho1)
#need to flip rotate 3 on yaxis
#extent(rotate3) <- extent(c(min(longitude-360),max(longitude-360),min(latitude),max(latitude)))
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am1 <- rotate3

##############################################################
###################################################################################################
# NOW AM2

all <- as.data.frame(cbind(jet$YEAR,jet[,11])) #indices_names[11]
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= 116 & longitude >=74)]
longitude <- longitude[which(longitude <= 116 & longitude >=74)]

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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am2 <- rotate3
###################################################################################################
# AM3
indices_names[12]
all <- as.data.frame(cbind(jet$YEAR,jet[,12])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= 152 & longitude >=116)]
longitude <- longitude[which(longitude <= 152 & longitude >=116)]

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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am3 <- rotate3
###################################################################################################
# AM4 #part 1 
q <- 13
lonmin1 <- 152
lonmax1 <- 180

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am4 <- rotate3
############################################################################
###################################################################################################
# AM4 # part 2 
q <- 13
lonmin1 <- -180
lonmax1 <- -150

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am4b <- rotate3

###################################################################################################
# AM5
q <- 14
lonmin1 <- -150
lonmax1 <- -120

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

am5 <- rotate3

###################################################################################################

# AM6 
q <- 15
lonmin1 <- -120
lonmax1 <- -94

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
#rotate3 <- raster(rho1) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

am6 <- rotate3
###################################################################################################
# AM7
q <- 16
lonmin1 <- -94
lonmax1 <- -56

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

am7 <- rotate3

###################################################################################################
# AM8 
q <- 17
lonmin1 <- -10
lonmax1 <- 8

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

pdsi2 <- pdsi1[which(1901:2017 <=2012 & 1901:2017 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
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
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))


am8 <- rotate3
##########################################################################################################################
t1 <- merge(am1,am2,am3,tolerance = 0.5)
t <- merge(am4,am4b,am5,am6,am7,am8, tolerance = 0.5) 

bg <- brewer.pal(n = 7, name = 'BrBG')

par(mai=c(0,0,0,0))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 


for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg1)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(24:74,rep(S[i],length(24:74)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg2)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(74:116,rep(S[i],length(116:74)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg3)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(116:152,rep(S[i],length(116:152)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg4)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(152:180,rep(S[i],length(152:180)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg4)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-180:-150,rep(S[i],length(-180:-150)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg5)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-150:-120,rep(S[i],length(-150:-120)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg6)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-120:-94,rep(S[i],length(-120:-94)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg7)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-94:-56,rep(S[i],length(-94:-56)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg8)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-10:8,rep(S[i],length(-10:8)), col = "grey")
}


abline(v = 24)
abline(v = 74)
abline(v = 116)
abline(v = 152)
abline(v = -150)
abline(v = -120)
abline(v = -94)
abline(v = -56)
abline(v = -10) 
abline(v = 8) 

map.axes()

plot(t1,main="",col =bg,breaks = c(-2.5,-1.5,-1,-0.5,0.5,1,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t,main="",col =bg,breaks = c(-2.5,-1.5,-1,-0.5,0.5,1,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)


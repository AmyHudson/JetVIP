# Jet Stream Hypothesis Figure GPH

setwd("/Volumes/AOP-NEON1.4/VIP/")

#Load Jet stream


# Read in AM Jet stream 
jet <- read.table("JetIndices.txt", header = TRUE)
jet <- jet[jet$YEAR>=1981,]
indices_names <- colnames(jet)
jet[,indices_names == "AM_Reg1"] #jet[,10]

library(maps)
library(mapdata)
library(maptools)  #for shapefiles
library(scales)  #for transparency
library(ggplot2)
library(devtools)
library(dplyr)
library(stringr)
library(ggmap) #need to cite ggmap

latmax <- 70
latmin <- 20
lonmax <- 180
lonmin <- -180


##################################################################################

# Now adding composite relationship with temperature field for subsetted region.

setwd("/Volumes/AOP-NEON1.4/monarch_transport/")

library(ncdf4)
library(fields)
library(Hmisc)
library(mapdata)
library(raster)
# compare bighorns instrumental june july august in bighorns grid cell with temperature in north america

#read in temperature field
data_time <- 1901:2014
#1901:2016 monthly
year <- matrix(rep(1851:2014,12),nrow = length(1851:2014),ncol = 12)
year <- as.vector(t(year))
month <- rep(1:12,length(1851:2014))#repeat 1:12 116 times
#repeat 1901:2016 12 times and make an array for the year vector
#crop to 1994 (jan or feb)
time <- as.data.frame(cbind(year,month))


# Monthly values for 1851/01 to 2014/12
b <- nc_open("hgt.mon.mean.nc") #91  24  180 1968

xdim <- round(b$dim[[3]]$vals, digits = 5)# lon 0 to 358
ydim <- round(b$dim[[1]]$vals, digits = 5) # lat 90 to -90
ldim <- b$dim[[2]]$vals #level
zdim <- b$dim[[4]]$vals # time

#keeping all longitdues, subsetting latitudes 20 to 70
xs <- which(xdim == 0) 
ys <-  which(ydim == 70) # 20.75
zs <- which(time$year == 1901)[1]
#zs <- which(zdim == 34348) # #know this from looking up the time variable 1117

hgt <- b$var[[2]] # lon 0 358 lat 90 to -90 level 11; 500hpa time 1968
hgt <- b$var$hgt # 180   91   24 1968 #x y l z

Clim_var.nc<- ncvar_get(b, hgt,start = c(xs,ys,1,zs), count = c(length(xdim)-xs+1,(which(ydim == 20)-which(ydim == 70)+1),length(ldim),length(zdim)-zs+1))    

dim(Clim_var.nc)[1]
dim(Clim_var.nc)[2]
dim(Clim_var.nc)[3]
dim(Clim_var.nc)[4]

Clim_var.nc <- Clim_var.nc[,,11,] #11 subsets 500hpa; 15 is 300hpa

time <- subset(time,year >=1901)
month <- time$month
year <- unique(time$year)

longitude <- xdim[c(xs:(length(xdim)))]
latitude <- ydim[c(ys:(ys+(which(ydim == 20)-which(ydim == 70))))]


mon1 <- Clim_var.nc[,,month == 4] #not sig april may
mon2 <- Clim_var.nc[,,month == 5]
mon <- (mon1 + mon2) / 2

tmp1 <- aperm(mon,c(3,2,1)) #reorder with time variable in front, lat, lon
#time is now in years
#mon1 <- aperm(mon,c(3,2,1)) #reorder with time variable in front, lat, lon
#dim(mon1)<- c(length(data_time), length(latitude)*length(longitude)) #reshape 2D w time, lat*lon
#mon1df <- as.data.frame(mon1)

library(dplyr)
year <- 1981:2012
###############################################################################
# Grab 10 south jet streams Region 1: 24 to 74

all <- as.data.frame(cbind(jet$YEAR,jet[,10]))
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 74 & longitude >=24)]
#crop to 1981:2012
#crop to 24 to 74 AM_1
longitude <- longitude[which(longitude <= 74 & longitude >=24)]
#latitude <- latitude[which(latitude <= 69.25 & latitude >=24.25)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) #reshape 2D w time, lat*lon
tmp2df <- as.data.frame(tmp2)



#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

my.t.test.p.value <- function(...) { 
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T) 
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

am1 <- rotate3
plot(am1)
####################################################################################################
# NOW AM2

all <- as.data.frame(cbind(jet$YEAR,jet[,11])) #indices_names[11]
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 116 & longitude >=74)]
longitude <- longitude[which(longitude <= 116 & longitude >=74)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], na.omit = T)
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
plot(am2)
###################################################################################################
# AM3
indices_names[12]
all <- as.data.frame(cbind(jet$YEAR,jet[,12])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 152 & longitude >=116)]
longitude <- longitude[which(longitude <= 152 & longitude >=116)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
plot(am3)
###################################################################################################
# AM4 #part 1 
q <- 13
lonmin1 <- 152
lonmax1 <- 180

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
plot(am4)
############################################################################
###################################################################################################
# AM4 # part 2 
q <- 13
lonmin1 <- 360-180
lonmax1 <- 360-150

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
extent(rotate3) <- extent(c(min(longitude)-360,max(longitude)-360,min(latitude),max(latitude)))

am4b <- rotate3
plot(am4b)
###################################################################################################
# AM5
q <- 14
lonmin1 <- 360-150
lonmax1 <- 360-120

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
extent(rotate3) <- extent(c(min(longitude)-360,max(longitude)-360,min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

am5 <- rotate3
plot(am5)
###################################################################################################
# AM6 
q <- 15
lonmin1 <- 360-120
lonmax1 <- 360-94

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
extent(rotate3) <- extent(c(min(longitude)-360,max(longitude)-360,min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

am6 <- rotate3
plot(am6)
###################################################################################################
# AM7
q <- 16
lonmin1 <- 360-94
lonmax1 <- 360-56

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
extent(rotate3) <- extent(c(min(longitude)-360,max(longitude)-360,min(latitude),max(latitude)))

am7 <- rotate3
plot(am7)
################################################################################################### LEFT OFF HERE_ split longitude

# AM8 
q <- 17
lonmin1 <- 360-10
lonmax1 <- 358

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
extent(rotate3) <- extent(c(min(longitude)-360,max(longitude)-360,min(latitude),max(latitude)))


am8 <- rotate3
plot(am8)
##########################################################################################################################
# AM8b 
q <- 17
lonmin1 <- 0
lonmax1 <- 8

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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


am8b <- rotate3
plot(am8b)
##########################################################################################################################


t <- merge(am5,am6,am7, tolerance = 0.5) 
# 

rwb <- (brewer.pal(n = 11, name = 'RdBu'))

png("/Volumes/AOP-NEON1.4/monarch_transport/gphAM_hi2.png",13,8,
    units = "in",res = 600, pointsize=20, family= "helvetica") 

par(mai=c(1,1,1,1))
map("world", xlim=c(-130,-50),ylim=c(20,70), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F)
map("world", xlim=c(-130,-50),ylim=c(20,70), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F)
#map.axes()


#plot(am6_diff_all,col = rwb, add = T, legend = F)#,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), add = T, legend = T)
rb <- rev(brewer.pal(n = 6, name = 'RdBu'))

plot(t,main="",col =rb,   xlab="",ylab="",horizontal = F,add = T, legend = T, breaks = c(-45,-30,-15,0,15,30,45))
#map("world", xlim=c(-130,-50),ylim=c(36,55), fill = F,  col="black", add = T)#, wrap=c(-180,180),add = F)
plot(t,main="",col =rb,breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)
dev.off()



#abline(v = -150)
abline(v = -120)
abline(v = -94)
abline(v = -56)



par(mai=c(0,0,0,0))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
#map("world", xlim=c(-180,180),ylim=c(20,70), fill = F,  col="black", add = F)#, wrap=c(-180,180),add = F)
#map("world", xlim=c(-180,180),ylim=c(20,70), fill = F,  col="black", add = T)#, wrap=c(-180,180),add = F)

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

t1 <- merge(am1,am2,am3,tolerance = 0.5)
t <- merge(am4,am4b,am5,am6,am7, tolerance = 0.5) 

library(RColorBrewer)
rb <- rev(brewer.pal(n = 6, name = 'RdBu'))

plot(t,main="",col =rb,   xlab="",ylab="",horizontal = F,add = T, legend = F, breaks = c(-45,-30,-15,0,15,30,45))
plot(t1,main="",col =rb,   xlab="",ylab="",horizontal = F,add = T, legend = F,breaks = c(-45,-30,-15,0,15,30,45)) #breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5),

t2 <- merge(am8,am8b, tolerance = 0.5) 


#t <- extend(t,e)
# otherwise I get a warning message 'number of items to replace is not a multiple of replacement length'
#t2 <- merge(am1,am2,am3,am4,am4b,am5,am6,am7,am8,tolerance = 0.7)

# Plot all 
#par(mai=c(0,0,0,0))  
#
#rb <- rev(brewer.pal(n = 5, name = 'RdBu'))#the color bar is still off here.

#map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(t1,main="",col =rb,breaks = c(-1.25,-1,-0.75,-0.5,0.5,0.75,1,1.25),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t,main="",col =rb,breaks = c(-1.25,-1,-0.75,-0.5,0.5,0.75,1,1.25),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t,main="",col =rb,breaks = c(-1.25,-1,-0.75,-0.5,0.5,0.75,1,1.25),   xlab="",ylab="",horizontal = F,add = T, legend = F)


plot(t1,main="",col =rb,breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t,main="",col =rb,breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)

# 
# 
# plot(am1,main="",col =rb,breaks = c(-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1),   xlab="",ylab="",horizontal = F,add = T, legend = F)
# plot(am2,main="",col =rb,breaks = c(-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1),   xlab="",ylab="",horizontal = F,add = T, legend = F)
# plot(am3,main="",col =rb,breaks = c(-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1),   xlab="",ylab="",horizontal = F,add = T, legend = F)
# plot(am4,main="",col =rb,breaks = c(-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1),   xlab="",ylab="",horizontal = F,add = T, legend = F)
# plot(am4b,main="",col =rb,breaks = c(-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1),   xlab="",ylab="",horizontal = F,add = T, legend = F)
# plot(am5,main="",col =rb,breaks = c(-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1),   xlab="",ylab="",horizontal = F,add = T, legend = F)
# plot(am6,main="",col =rb,breaks = c(-1.25,-1,-0.75,-0.25,0.25,0.75,1,1.25),   xlab="",ylab="",horizontal = F,add = T, legend = T)
# plot(am7,main="",col =rb,breaks = c(-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1),   xlab="",ylab="",horizontal = F,add = T, legend = F)
# plot(am8,main="",col =rb,breaks = c(-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1),   xlab="",ylab="",horizontal = F,add = T, legend = F)


# The plots are only correct spatially for the first two 'add = T' rasters... no matter which 2 rasters... perhaps because this is too much information to be stored? 
# if so, I will merge all the plots into one raster.

#extend function?
# r <- raster(xmn=-150, xmx=-120, ymx=60, ymn=30, ncol=36, nrow=18)
# r[] <- 1:ncell(r)
# e <- extent(-180, 180, 24, 70)
# re <- extend(r, e)
# plot(r)
# plot(re)


# 
# rasterOptions(tolerance = 0.1)
# .Machine$double.eps <- 0.000000001
# 
# library(raster)
# dem_n1<-getData("SRTM",lon=26,lat=43)
# dem_n2<-getData("SRTM",lon=31,lat=43)
# demALL<-merge(dem_n1,dem_n2)

# t <- merge(am1,am2,am3,am4,am4b,am5,am6,am7,am8, tolerance = 0.5)
# plot(t)

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

# Downloaded Countries outline in the form of ne_50m_admin_0_countries.zip from: https://www.naturalearthdata.com/downloads/50m-cultural-vectors/50m-admin-0-countries-2/
#7/11/2018
require(utils)
require(RNetCDF)
#require(rasterVis)
library(rasterVis)
library(rgdal)

#cntry <- readOGR('./ne_50m_admin_0_countries.shp',stringsAsFactors = TRUE)


require(OceanView) 
library(plot3D)
library(lattice)

rb <- rev(brewer.pal(n = 7, name = 'RdBu'))#the color bar is still off here.

# correlations 
# setEPS()
# postscript("bhjja.eps")
# plot(rotate3,main="",col =rb,breaks = c(-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1),   xlab="",ylab="",horizontal = F)
# #map("world",add=T,lwd=2)
# map("worldHires","Canada", xlim=c(-150,max(longitude)),ylim=c(min(latitude),max(latitude)), fill=F,lwd = 2, add = T)  #plot the region of Canada
# map("worldHires","USA", xlim=c(-150,max(longitude)),ylim=c(min(latitude),max(latitude)), fill=F, add = T,lwd = 2)  #plot the region of US
# map("worldHires","Mexico", xlim=c(-150,max(longitude)),ylim=c(min(latitude),max(latitude)), fill=F, add = T,lwd = 2)#plot the region of Mexico
# points(-107.7, 44.5,cex = 1.5, pch=19, col="black")  
# mtext("A", side=3, adj=0, line=0.5, cex=2, font=2); 
# dev.off()

#####################################################################

# Grab 10 NORTH jet streams Region 1: 24 to 74
longitude <- xdim[c(xs:(length(xdim)))]


all <- as.data.frame(cbind(jet$YEAR,jet[,10]))
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 74 & longitude >=24)]
#crop to 1981:2012
#crop to 24 to 74 AM_1
longitude <- longitude[which(longitude <= 74 & longitude >=24)]
#latitude <- latitude[which(latitude <= 69.25 & latitude >=24.25)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) #reshape 2D w time, lat*lon
tmp2df <- as.data.frame(tmp2)



m1 <- tmp2df %>% filter(year %in% as.character(allN))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

my.t.test.p.value <- function(...) { 
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T) #greater for jet poleward
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
plot(am1)
####################################################################################################
# NOW AM2

all <- as.data.frame(cbind(jet$YEAR,jet[,11])) #indices_names[11]
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 116 & longitude >=74)]
longitude <- longitude[which(longitude <= 116 & longitude >=74)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
plot(am2)
###################################################################################################
# AM3
indices_names[12]
all <- as.data.frame(cbind(jet$YEAR,jet[,12])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= 152 & longitude >=116)]
longitude <- longitude[which(longitude <= 152 & longitude >=116)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
plot(am3)
###################################################################################################
# AM4 #part 1 
q <- 13
lonmin1 <- 152
lonmax1 <- 180

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
plot(am4)
############################################################################
###################################################################################################
# AM4 # part 2 
q <- 13
lonmin1 <- 360-180
lonmax1 <- 360-150

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
extent(rotate3) <- extent(c(min(longitude)-360,max(longitude)-360,min(latitude),max(latitude)))

am4b <- rotate3
plot(am4b)
###################################################################################################
# AM5
q <- 14
lonmin1 <- 360-150
lonmax1 <- 360-120

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
extent(rotate3) <- extent(c(min(longitude)-360,max(longitude)-360,min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

am5 <- rotate3
plot(am5)
###################################################################################################
# AM6 
q <- 15
lonmin1 <- 360-120
lonmax1 <- 360-94

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
extent(rotate3) <- extent(c(min(longitude)-360,max(longitude)-360,min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

am6 <- rotate3
plot(am6)
###################################################################################################
# AM7
q <- 16
lonmin1 <- 360-94
lonmax1 <- 360-56

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
extent(rotate3) <- extent(c(min(longitude)-360,max(longitude)-360,min(latitude),max(latitude)))

am7 <- rotate3
plot(am7)
###################################################################################################
# AM8 
q <- 17
lonmin1 <- 360-10
lonmax1 <- 358

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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
extent(rotate3) <- extent(c(min(longitude)-360,max(longitude)-360,min(latitude),max(latitude)))


am8 <- rotate3
plot(am8)
##########################################################################################################################
###############################################################
# AM8b 
q <- 17
lonmin1 <- 0
lonmax1 <- 8

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(length(xdim)))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))
#m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i],na.omit = T)
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


am8b <- rotate3
plot(am8b)
################
par(mai=c(0,0,0,0))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

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

rb <- rev(brewer.pal(n = 6, name = 'RdBu'))

t1 <- merge(am1,am2,am3,tolerance = 0.5) 
t <- merge(am4,am4b,am5,am6,am7,tolerance = 0.5) 
t <- merge(am5,am6,am7,tolerance = 0.5) 

plot(t,main="",col =rb,breaks = c(-45,-30,-15,0,15,30,45),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t1,main="",col =rb,breaks = c(-45,-30,-15,0,15,30,45),   xlab="",ylab="",horizontal = F,add = T, legend = F)

#JetHypothesisFig_JF_temp.R


#Load Jet stream


# Read in Jet stream 
jet <- read.table("JetIndices.txt", header = TRUE)
jet <- jet[jet$YEAR>=1981,]
indices_names <- colnames(jet)
which(indices_names == "JF_Reg1") #jet[,19]

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

setwd("/Volumes/AOP-NEON1.4/monarch_transport/")

library(ncdf4)
library(fields)
library(Hmisc)
library(mapdata)

#read in temperature field
data_time <- 1901:2016
#1901:2016 monthly
year <- matrix(rep(1901:2016,12),nrow = length(1901:2016),ncol = 12)
year <- as.vector(t(year))
month <- rep(1:12,116)#repeat 1:12 116 times
#repeat 1901:2016 12 times and make an array for the year vector
#crop to 1994 (jan or feb)
time <- cbind(year,month)

a <- nc_open("cru_ts4.01.1901.2016.tmp.dat.nc")

xdim <- round(a$dim[[1]]$vals, digits = 5)# lon
ydim <- round(a$dim[[2]]$vals, digits = 5) # lat
zdim <- a$dim[[3]]$vals # time

xs <- which(xdim == -179.75)
ys <-  which(ydim == 24.25) #allows us to crop to midlatitudes
#zs <- which(zdim == 34348) # #know this from looking up the time variable 1117

tmp <- a$var[[1]]

ptm <- proc.time()   
#Clim_var.nc<- ncvar_get(a, tmp,start = c(xs,ys,zs), count = c(114, 92,276))  
Clim_var.nc<- ncvar_get(a, tmp,start = c(xs,ys,1), count = c((which(xdim == 179.75)-which(xdim == -179.75)+1), (which(ydim == 69.25)-which(ydim == 24.25)+1),1392))    
proc.time() - ptm
dim(Clim_var.nc)[1]
dim(Clim_var.nc)[2]
dim(Clim_var.nc)[3]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]
latitude <- ydim[c(ys:(ys+(which(ydim == 69.25)-which(ydim == 24.25))))]

mon1 <- Clim_var.nc[,,month == 1]
mon2 <- Clim_var.nc[,,month == 2]
mon <- (mon1 + mon2) / 2

tmp1 <- aperm(mon,c(3,2,1)) #reorder with time variable in front, lat, lon
#time is now in years
###########################################################################################
# Grab 10 north/south jet streams 
#JF1
q <- 2
lonmin1 <- 18
lonmax1 <- 40

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) #reshape 2D w time, lat*lon
tmp2df <- as.data.frame(tmp2)

year <- 1981:2012

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

my.t.test.p.value <- function(...) { 
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "less",na.omit = T) #greater for jet poleward
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
rotate3 <- raster(rho1[nrow(rho1):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j1 <- rotate3

####################################################################################################
# JF2
q <- 3
lonmin1 <- 40
lonmax1 <- 80

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "less",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j2 <- rotate3
map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j2,add = T)
###################################################################################################
# JF3
q <- 4
lonmin1 <- 80
lonmax1 <- 146

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "less",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j3 <- rotate3
map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 


plot(j3,add = T)
###################################################################################################
# JF4 
q <- 5
lonmin1 <- 146
lonmax1 <- 172

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "less",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j4 <- rotate3
############################################################################
###################################################################################################
# JF5
q <- 6
lonmin1 <- 172
lonmax1 <- 180

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "less",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j5 <- rotate3

###################################################################################################
# JF5
q <- 6
lonmin1 <- -180
lonmax1 <- -100

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "less",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

j5b <- rotate3

###################################################################################################
# JF6  
q <- 7
lonmin1 <- -100
lonmax1 <- -72

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "less",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

j6 <- rotate3

###################################################################################################
# JF7
q <- 8
lonmin1 <- -72
lonmax1 <- -28

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "less",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j7 <- rotate3
###################################################################################################
# JF8
q <- 9
lonmin1 <- -28
lonmax1 <- 8

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

#m1 <- tmp2df %>% filter(year %in% as.character(allN))
m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "less",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j8 <- rotate3
###################################################################################################

##########################################################################################################################
t1 <- merge(j2,j1,tolerance = 0.5)
map("world", xlim=c(10,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
#map("world", xlim=c(20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
plot(t1,add = T,legend = F)
plot(j3,add = T,legend = F)
#plot(j4,add = T,legend = F)
#plot(j5,add = T,legend = F)

t2 <- merge(j7,j8,tolerance = 0.5)

#nothing in j6
#t <- merge(j5b,j6, tolerance = 0.5) 
map("world", xlim=c(-180,8),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
#t2 <- merge(j4,j5, tolerance = 0.5) 
plot(t2,add = T, legend = F)
plot(j6,add = T, legend = F)
#plot(j5b,add = T, legend = F)

#plot(j5b,add = T, legend = F)

# Plot all 

par(mai=c(0,0,0,1))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 


for (i in 1:length(jet$JF_Reg1)){
  lines(18:40,rep(jet$JF_Reg1[i],length(18:40)), col = "grey")
}
for (i in 1:length(jet$JF_Reg2)){
  lines(40:80,rep(jet$JF_Reg2[i],length(40:80)), col = "grey")
}
for (i in 1:length(jet$JF_Reg3)){
  lines(80:146,rep(jet$JF_Reg3[i],length(80:146)), col = "grey")
}
for (i in 1:length(jet$JF_Reg4)){
  lines(146:172,rep(jet$JF_Reg3[i],length(146:172)), col = "grey")
}
for (i in 1:length(jet$JF_Reg5)){
  lines(172:180,rep(jet$JF_Reg4[i],length(172:180)), col = "grey")
}
for (i in 1:length(jet$JF_Reg5)){
  lines(-180:-100,rep(jet$JF_Reg5[i],length(-180:-100)), col = "grey")
}
for (i in 1:length(jet$JF_Reg6)){
  lines(-100:-72,rep(jet$JF_Reg6[i],length(-100:-72)), col = "grey")
}
for (i in 1:length(jet$JF_Reg7)){
  lines(-72:-28,rep(jet$JF_Reg7[i],length(-72:-28)), col = "grey")
}
for (i in 1:length(jet$JF_Reg8)){
  lines(-28:8,rep(jet$JF_Reg7[i],length(-28:8)), col = "grey")
}

abline(v = 18)
abline(v = 40)
abline(v = 80)
abline(v = 146)
abline(v = 172)
abline(v = -100)
abline(v = -72)
abline(v = -28)
abline(v = 8)

map.axes() 


rb <- rev(brewer.pal(n = 7, name = 'RdBu'))
rb <- c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7","#F7F7F7", "#FDDBC7", "#EF8A62","#B2182B")


par(mai=c(0,0,0,0))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$ON_Reg1)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(18:40,rep(S[i],length(18:40)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg2)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(40:80,rep(S[i],length(40:80)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg3)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(80:146,rep(S[i],length(80:146)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg4)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(146:172,rep(S[i],length(146:172)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg5)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(172:180,rep(S[i],length(172:180)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg5)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-180:-100,rep(S[i],length(-180:-100)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg6)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-100:-72,rep(S[i],length(-100:-72)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg7)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-72:-28,rep(S[i],length(-72:-28)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg8)) #
  S <- all[order(all$V2),]$V2[1:10]
  lines(-28:8,rep(S[i],length(-28:8)), col = "grey")
}




abline(v = 18)
abline(v = 40)
abline(v = 80)
abline(v = 146)
abline(v = 172)
abline(v = -100)
abline(v = -72)
abline(v = -28)
abline(v = 8)

map.axes()

plot(t1,main="",col =rb,breaks = c(-2,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j3,main="",col =rb,breaks = c(-2,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)

#plot(j4,main="",col =rb,breaks = c(-2,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)
#plot(j5,main="",col =rb,breaks = c(-2,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)

plot(t2,main="",col =rb,breaks = c(-2,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j6,main="",col =rb,breaks = c(-2,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,2),   xlab="",ylab="",horizontal = F,add = T, legend = F)

##############################################################################################################################

## NOW POLEWARD!

#JF1
q <- 2
lonmin1 <- 18
lonmax1 <- 40

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) #reshape 2D w time, lat*lon
tmp2df <- as.data.frame(tmp2)

year <- 1981:2012

m1 <- tmp2df %>% filter(year %in% as.character(allN))
#m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

my.t.test.p.value <- function(...) { 
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "greater",na.omit = T) #greater for jet poleward
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
rotate3 <- raster(rho1[nrow(rho1):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j1 <- rotate3

####################################################################################################
# JF2
q <- 3
lonmin1 <- 40
lonmax1 <- 80

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))
#m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "greater",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j2 <- rotate3
map("world", xlim=c(54,86),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j2,add = T)
###################################################################################################
# JF3
q <- 4
lonmin1 <- 80
lonmax1 <- 146

indices_names[q]

all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >= lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >= lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))
#m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "greater",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j3 <- rotate3
#map("world", xlim=c(86,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

#plot(j3,add = T)
###################################################################################################
# JF4 
q <- 5
lonmin1 <- 146
lonmax1 <- 172

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))
#m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "greater",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j4 <- rotate3
############################################################################
###################################################################################################
# JF5
q <- 6
lonmin1 <- 172
lonmax1 <- 180

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))
#m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "greater",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j5 <- rotate3

###################################################################################################
# JF5
q <- 6
lonmin1 <- -180
lonmax1 <- -100

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))
#m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "greater",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

j5b <- rotate3

###################################################################################################
# JF6 
q <- 7
lonmin1 <- -100
lonmax1 <- -72

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))
#m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "greater",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))
#extent(rotate3) <- extent(c(-180,180,24,70))

j6 <- rotate3

###################################################################################################
# JF7
q <- 8
lonmin1 <- -72
lonmax1 <- -28

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))
#m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "greater",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j7 <- rotate3
###################################################################################################
# JF8
q <- 9
lonmin1 <- -28
lonmax1 <- 8

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

longitude <- xdim[c(xs:(xs+(which(xdim == 179.75)-which(xdim == -179.75))))]

tmp2 <- tmp1[which(1901:2016 <=2012 & 1901:2016 >=1981),,which(longitude <= lonmax1 & longitude >=lonmin1)]
longitude <- longitude[which(longitude <= lonmax1 & longitude >=lonmin1)]

dim(tmp2)<- c(length(1981:2012),length(latitude)*length(longitude)) 
tmp2df <- as.data.frame(tmp2)

m1 <- tmp2df %>% filter(year %in% as.character(allN))
#m1 <- tmp2df %>% filter(year %in% as.character(allS))

diff <- colMeans(m1)-colMeans(tmp2df) #a positive means warmer

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(diff)))#NA[length(mon1)]

for (i in 1:length(diff)){
  y <- my.t.test.p.value(m1[,i],tmp2df[,i], alternative = "greater",na.omit = T)
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
rotate3 <- raster(rho1[nrow(rho1):1,]) 
extent(rotate3) <- extent(c(min(longitude),max(longitude),min(latitude),max(latitude)))

j8 <- rotate3

###################################################################################################

##########################################################################################################################
t1 <- merge(j1,j2,tolerance = 0.5)
#map("world", xlim=c(20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(-20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
plot(t1,add = T,legend = F)
plot(j3,add = T,legend = F)

#t2 <- merge(j4,j5,tolerance = 0.5)

#nothing in j6
t <- merge(j6,j8, tolerance = 0.5) 
map("world", xlim=c(-180,8),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
#t2 <- merge(j4,j5, tolerance = 0.5) 
plot(j5b,add = T, legend = F)
plot(t,add = T, legend = F)
#plot(j8,add = T, legend = F)

# Plot all 


rb <- rev(brewer.pal(n = 7, name = 'RdBu'))
rb <- c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7","#F7F7F7", "#FDDBC7", "#EF8A62","#B2182B")


par(mai=c(0,0,0,1))  
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(lonmin,lonmax),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg1)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(18:40,rep(S[i],length(18:40)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg2)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(40:80,rep(S[i],length(40:80)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg3)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(80:146,rep(S[i],length(80:146)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg4)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(146:172,rep(S[i],length(146:172)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg5)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(172:180,rep(S[i],length(172:180)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg5)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-180:-100,rep(S[i],length(-180:-100)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg6)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-100:-72,rep(S[i],length(-100:-72)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg7)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-72:-28,rep(S[i],length(-72:-28)), col = "grey")
}

for (i in 1:10){
  all <- as.data.frame(cbind(jet$YEAR,jet$JF_Reg8)) #
  S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
  lines(-28:8,rep(S[i],length(-28:8)), col = "grey")
}

abline(v = 18)
abline(v = 40)
abline(v = 80)
abline(v = 146)
abline(v = 172)
abline(v = -100)
abline(v = -72)
abline(v = -28)
abline(v = 8)

map.axes()

plot(t1,main="",col =rb,breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j3,main="",col =rb,breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),   xlab="",ylab="",horizontal = F,add = T, legend = F)

plot(j5b,main="",col =rb,breaks = c(-1.5,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t,main="",col =rb,breaks = c(-1.5,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)

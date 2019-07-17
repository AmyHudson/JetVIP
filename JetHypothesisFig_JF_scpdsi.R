#JetHypothesisFig_JF_scpdsi.R

# Jet Stream Hypothesis Figure: I envision April May Jet stream as faint lines over the northern hemisphere background, with extreme poleward (10) jets bolded and extreme equatorward bolded. 

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

mon1 <- Clim_var.nc[,,month == 1] 
mon2 <- Clim_var.nc[,,month == 2]
mon <- (mon1 + mon2) / 2

pdsi1 <- aperm(mon,c(3,2,1)) #reorder with time variable in front, lat, lon
#time is now in years

###########################################################################################
# Grab 10 north/south jet streams ON Region 1
#JF1
q <- 2
lonmin1 <- 18
lonmax1 <- 40

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
# JF2
q <- 3
lonmin1 <- 40
lonmax1 <- 80

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
# JF3
q <- 4
lonmin1 <- 80
lonmax1 <- 146

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
# JF4 
q <- 5
lonmin1 <- 146
lonmax1 <- 172

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
# JF5
q <- 6
lonmin1 <- 172
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
# JF5
q <- 6
lonmin1 <- -180
lonmax1 <- -100

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
# JF6 
q <- 7
lonmin1 <- -100
lonmax1 <- -72


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
# JF7
q <- 8
lonmin1 <- -72
lonmax1 <- -28

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
# JF8
q <- 9
lonmin1 <- -28
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

j8 <- rotate3

map("world", xlim=c(lonmin1,lonmax1),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 

plot(j8,add = T)
##########################################################################################################################
################################################################################################################
t1 <- merge(j1,j2,tolerance = 0.5)
#t3 <- merge(j3,j4,tolerance = 0.5)
#map("world", xlim=c(20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
#plot(t3,add = T,legend = F)
plot(t1,add = T,legend = F)
plot(j3,add = T,legend = F)
plot(j4,add = T,legend = F)


#nothing in j6
#t <- merge(j6, tolerance = 0.5) 
map("world", xlim=c(-180,8),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
t2 <- merge(j7,j8, tolerance = 0.5) 

plot(t2,add = T, legend = F)
plot(j5b,add = T, legend = F)
plot(j6,add = T, legend = F)


# Plot all 

bg <- brewer.pal(n = 7, name = 'BrBG')
bg <- brewer.pal(n = 5, name = 'BrBG')

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

plot(t1,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j3,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)

plot(j4,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)

plot(j5b,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j6,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(t2,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)


########################################################################################################################

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
# JF2
q <- 3
lonmin1 <- 40
lonmax1 <- 80

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
# JF4 
q <- 5
lonmin1 <- 146
lonmax1 <- 172

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
################################################################# JF5
q <- 6
lonmin1 <- 172
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
# JF5
q <- 6
lonmin1 <- -180
lonmax1 <- -100

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

################################################################# JF6 
q <- 7
lonmin1 <- -100
lonmax1 <- -72

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
# JF7
q <- 8
lonmin1 <- -72
lonmax1 <- -28


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
# JF8
q <- 9
lonmin1 <- -28
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

j8 <- rotate3
###################################################################################################
t1 <- merge(j1,j2,tolerance = 0.5)
#map("world", xlim=c(20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
map("world", xlim=c(-20,180),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
plot(t1,add = T,legend = T)
plot(j3,add = T,legend = T)
plot(j4,add = T,legend = T)
plot(j5,add = T,legend = T)

#nothing in j6
t <- merge(j7,j8, tolerance = 0.5) 
map("world", xlim=c(-180,8),ylim=c(latmin,latmax), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F) 
#t2 <- merge(j4,j5, tolerance = 0.5) 
plot(t,add = T, legend = F)
# Plot all 
plot(j5b,add = T, legend = T)
plot(j6,add = T, legend = T)


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


plot(t1,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j3,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j4,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j5,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)


plot(t,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j5b,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)
plot(j6,main="",col =bg,breaks = c(-2.5,-1.5,-0.5,0.5,1.5,2.5),   xlab="",ylab="",horizontal = F,add = T, legend = F)

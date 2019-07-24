
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
library(raster)

library(ncdf4)
library(fields)
library(Hmisc)
library(mapdata)

my.t.test.p.value <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
my.t.test.statistic <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$statistic)
}


library(R.matlab)

year <- 1981:2014
lat <- seq(20,70,0.05)


####################################################################################################
#AM_1 
q <- 10
lonmin1 <- 24
lonmax1 <- 74
a <- readMat("SOS1_AM_R3_20N70N_1.mat")
lon <- seq(lonmin1,lonmax1,0.05)

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

a2 <- a$sos3[which(year <=2012 & year >=1981),,]
dim(a2)<- c(length(1981:2012),length(lat)*length(lon)) 
a2df <- as.data.frame(a2)

m1 <- a2df[year %in% as.character(allS),]
m2 <- m1
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am1_diff_south <- rotate3

m1 <- a2df[year %in% as.character(allN),]
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am1_diff_north <- rotate3

diff <- colMeans(m1)-colMeans(m2)
diff <- as.matrix(diff)
diff[diff>12] <- 12
diff[-12>diff] <- -12
dim(diff) <- c(length(lat),length(lon))
diff <- as.data.frame(diff)
colnames(diff) <- lon
rownames(diff) <- lat
diff <- as.matrix(diff)
rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
am1_diff_all <- rotate3

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))#NA[length(mon1)]

for (i in 1:length(m1)){
  y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T) 
  z <- my.t.test.statistic(m1[,i],m2[,i],na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- z
  }
}

#write.csv(rho1, file = "r1.csv")

rho2 <- as.matrix(rho1)
dim(rho2) <- c(length(lat),length(lon))
rho2 <- as.data.frame(rho2)
rho2[rho2>0] <- 1
rho2[0>rho2] <- -1
colnames(rho2) <- lon
rownames(rho2) <- lat
rho2 <- as.matrix(rho2)
rotate3 <- raster(rho2[nrow(rho2):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))

am1_ns_sig <- rotate3
write.csv(rho2, file = "am1_ns_sig.csv")
# ##############################################################
# ###################################################################################################
# # NOW AM2
q <- 11
lonmin1 <- 74
lonmax1 <- 116
a <- readMat("SOS1_AM_R3_20N70N_2.mat")
lon <- seq(lonmin1,lonmax1,0.05)

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

a2 <- a$sos3[which(year <=2012 & year >=1981),,]
dim(a2)<- c(length(1981:2012),length(lat)*length(lon))
a2df <- as.data.frame(a2)

m1 <- a2df[year %in% as.character(allS),]
m2 <- m1
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am2_diff_south <- rotate3
# 
m1 <- a2df[year %in% as.character(allN),]
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am2_diff_north <- rotate3
# 
diff <- colMeans(m1)-colMeans(m2)
diff <- as.matrix(diff)
diff[diff>12] <- 12
diff[-12>diff] <- -12
dim(diff) <- c(length(lat),length(lon))
diff <- as.data.frame(diff)
colnames(diff) <- lon
rownames(diff) <- lat
diff <- as.matrix(diff)
rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
am2_diff_all <- rotate3

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))#NA[length(mon1)]

for (i in 1:length(m1)){
  y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T) 
  z <- my.t.test.statistic(m1[,i],m2[,i],na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- z
  }
}

#write.csv(rho1, file = "r2.csv")

rho2 <- as.matrix(rho1)
dim(rho2) <- c(length(lat),length(lon))
rho2 <- as.data.frame(rho2)
rho2[rho2>0] <- 1
rho2[0>rho2] <- -1
colnames(rho2) <- lon
rownames(rho2) <- lat
rho2 <- as.matrix(rho2)
rotate3 <- raster(rho2[nrow(rho2):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))

am2_ns_sig <- rotate3
write.csv(rho2, file = "am2_ns_sig.csv")
# ###################################################################################################
# AM3

q <- 12
lonmin1 <- 116
lonmax1 <- 152
a <- readMat("SOS1_AM_R3_20N70N_3.mat")
lon <- seq(lonmin1,lonmax1,0.05)

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

a2 <- a$sos3[which(year <=2012 & year >=1981),,]
dim(a2)<- c(length(1981:2012),length(lat)*length(lon))
a2df <- as.data.frame(a2)

m1 <- a2df[year %in% as.character(allS),]
m2 <- m1
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am3_diff_south <- rotate3
# 
m1 <- a2df[year %in% as.character(allN),]
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am3_diff_north <- rotate3
# 
diff <- colMeans(m1)-colMeans(m2)
diff <- as.matrix(diff)
diff[diff>12] <- 12
diff[-12>diff] <- -12
dim(diff) <- c(length(lat),length(lon))
diff <- as.data.frame(diff)
colnames(diff) <- lon
rownames(diff) <- lat
diff <- as.matrix(diff)
rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
am3_diff_all <- rotate3
# 


rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))#NA[length(mon1)]

for (i in 1:length(m1)){
  y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T) 
  z <- my.t.test.statistic(m1[,i],m2[,i],na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- z
  }
}

#write.csv(rho1, file = "r3.csv")

rho2 <- as.matrix(rho1)
dim(rho2) <- c(length(lat),length(lon))
rho2 <- as.data.frame(rho2)
rho2[rho2>0] <- 1
rho2[0>rho2] <- -1
colnames(rho2) <- lon
rownames(rho2) <- lat
rho2 <- as.matrix(rho2)
rotate3 <- raster(rho2[nrow(rho2):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))

am3_ns_sig <- rotate3
write.csv(rho2, file = "am3_ns_sig.csv")
# ###################################################################################################
# AM4 #part 1
q <- 13
lonmin1 <- 152
lonmax1 <- 179.95

a <- readMat("SOS1_AM_R3_20N70N_4.mat")
lon <- seq(lonmin1,lonmax1,0.05)

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

a2 <- a$sos3[which(year <=2012 & year >=1981),,]
dim(a2)<- c(length(1981:2012),length(lat)*length(lon))
a2df <- as.data.frame(a2)

m1 <- a2df[year %in% as.character(allS),]
m2 <- m1
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am4_diff_south <- rotate3
# 
m1 <- a2df[year %in% as.character(allN),]
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am4_diff_north <- rotate3
# 
diff <- colMeans(m1)-colMeans(m2)
diff <- as.matrix(diff)
diff[diff>12] <- 12
diff[-12>diff] <- -12
dim(diff) <- c(length(lat),length(lon))
diff <- as.data.frame(diff)
colnames(diff) <- lon
rownames(diff) <- lat
diff <- as.matrix(diff)
rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
am4_diff_all <- rotate3

plot(am4_diff_all)

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))#NA[length(mon1)]

for (i in 1:length(m1)){
  y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T) 
  z <- my.t.test.statistic(m1[,i],m2[,i],na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- z
  }
}

#write.csv(rho1, file = "r4a.csv")

rho2 <- as.matrix(rho1)
dim(rho2) <- c(length(lat),length(lon))
rho2 <- as.data.frame(rho2)
rho2[rho2>0] <- 1
rho2[0>rho2] <- -1
colnames(rho2) <- lon
rownames(rho2) <- lat
rho2 <- as.matrix(rho2)
rotate3 <- raster(rho2[nrow(rho2):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))

am4a_ns_sig <- rotate3
plot(am4a_ns_sig, col = rwb)
plot(am4_diff_all, col = rwb)

write.csv(rho2, file = "am4a_ns_sig.csv")

# ############################################################################
# ###################################################################################################
# AM4 # part 2
q <- 13
lonmin1 <- -180
lonmax1 <- -150

a <- readMat("SOS1_AM_R3_20N70N_4b.mat")
lon <- seq(lonmin1,lonmax1,0.05)

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

a2 <- a$sos3[which(year <=2012 & year >=1981),,]
dim(a2)<- c(length(1981:2012),length(lat)*length(lon))
a2df <- as.data.frame(a2)

m1 <- a2df[year %in% as.character(allS),]
m2 <- m1
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am4b_diff_south <- rotate3
# 
m1 <- a2df[year %in% as.character(allN),]
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am4b_diff_north <- rotate3
# 
diff <- colMeans(m1)-colMeans(m2)
diff <- as.matrix(diff)
diff[diff>12] <- 12
diff[-12>diff] <- -12
dim(diff) <- c(length(lat),length(lon))
diff <- as.data.frame(diff)
colnames(diff) <- lon
rownames(diff) <- lat
diff <- as.matrix(diff)
rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
am4b_diff_all <- rotate3

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))#NA[length(mon1)]

for (i in 1:length(m1)){
  y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T) 
  z <- my.t.test.statistic(m1[,i],m2[,i],na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- z
  }
}

#write.csv(rho1, file = "r4b.csv")

rho2 <- as.matrix(rho1)
dim(rho2) <- c(length(lat),length(lon))
rho2 <- as.data.frame(rho2)
rho2[rho2>0] <- 1
rho2[0>rho2] <- -1
colnames(rho2) <- lon
rownames(rho2) <- lat
rho2 <- as.matrix(rho2)
rotate3 <- raster(rho2[nrow(rho2):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))

am4b_ns_sig <- rotate3
write.csv(rho2, file = "am4b_ns_sig.csv")
# ###################################################################################################
# AM5
q <- 14
lonmin1 <- -150
lonmax1 <- -120

a <- readMat("SOS1_AM_R3_20N70N_5.mat") #
lon <- seq(lonmin1,lonmax1,0.05)

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

a2 <- a$sos3[which(year <=2012 & year >=1981),,]
dim(a2)<- c(length(1981:2012),length(lat)*length(lon))
a2df <- as.data.frame(a2)

m1 <- a2df[year %in% as.character(allS),]
m2 <- m1
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am5_diff_south <- rotate3
# 
m1 <- a2df[year %in% as.character(allN),]
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am5_diff_north <- rotate3
# 
diff <- colMeans(m1)-colMeans(m2)
diff <- as.matrix(diff)
diff[diff>12] <- 12
diff[-12>diff] <- -12
dim(diff) <- c(length(lat),length(lon))
diff <- as.data.frame(diff)
colnames(diff) <- lon
rownames(diff) <- lat
diff <- as.matrix(diff)
rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
am5_diff_all <- rotate3

plot(am5_diff_all)
# rho2 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))
# 
# for (i in 1:length(m1)){
#   y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T)
#   if(y < 0.1 & is.na(y) == "FALSE"){
#     rho2[i] <- mean(m1[,i])-mean(m2[,i])
#   }
# }
# #maybe the median would be better
# 
# rho1 <- as.matrix(rho2)
# dim(rho1) <- c(length(lat),length(lon))
# rho1 <- as.data.frame(rho1)
# rho1[rho1>20] <- 20
# rho1[-20>rho1] <- -20
# colnames(rho1) <- lon
# rownames(rho1) <- lat
# rho1 <- as.matrix(rho1)
# rotate3 <- raster(rho1[nrow(rho1):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# 
# am5_ns_sig <- rotate3

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))#NA[length(mon1)]

for (i in 1:length(m1)){
  y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T) 
  z <- my.t.test.statistic(m1[,i],m2[,i],na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- z
  }
}

#write.csv(rho1, file = "r5.csv")

rho2 <- as.matrix(rho1)
dim(rho2) <- c(length(lat),length(lon))
rho2 <- as.data.frame(rho2)
rho2[rho2>0] <- 1
rho2[0>rho2] <- -1
colnames(rho2) <- lon
rownames(rho2) <- lat
rho2 <- as.matrix(rho2)
rotate3 <- raster(rho2[nrow(rho2):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))

am5_ns_sig <- rotate3
write.csv(rho2, file = "am5_ns_sig.csv")
# 
# ###################################################################################################
# 
# AM6
q <- 15
lonmin1 <- -120
lonmax1 <- -94
a <- readMat("SOS1_AM_R3_20N70N_6.mat")
lon <- seq(lonmin1,lonmax1,0.05)

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

a2 <- a$sos3[which(year <=2012 & year >=1981),,]
dim(a2)<- c(length(1981:2012),length(lat)*length(lon))
a2df <- as.data.frame(a2)
#a2df[a2df == "NaN"] = "NA" 

m1 <- a2df[year %in% as.character(allS),]
m2 <- m1
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am6_diff_south <- rotate3

m1 <- a2df[year %in% as.character(allN),]
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am6_diff_north <- rotate3

diff <- colMeans(m1)-colMeans(m2)
diff <- as.matrix(diff)
diff[diff>12] <- 12
diff[-12>diff] <- -12
dim(diff) <- c(length(lat),length(lon))
diff <- as.data.frame(diff)
colnames(diff) <- lon
rownames(diff) <- lat
diff <- as.matrix(diff)
rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
am6_diff_all <- rotate3

plot(am6_diff_all)
# rho2 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))
# 
# for (i in 1:length(m1)){
#   y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T)
#   if(y < 0.1 & is.na(y) == "FALSE"){
#     rho2[i] <- mean(m1[,i])-mean(m2[,i])
#   }
# }
#maybe the median would be better

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))#NA[length(mon1)]

for (i in 1:length(m1)){
  y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T) 
  z <- my.t.test.statistic(m1[,i],m2[,i],na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- z
  }
}

#write.csv(rho1, file = "r6.csv")

rho2 <- as.matrix(rho1)
dim(rho2) <- c(length(lat),length(lon))
rho2 <- as.data.frame(rho2)
rho2[rho2>0] <- 1
rho2[0>rho2] <- -1
colnames(rho2) <- lon
rownames(rho2) <- lat
rho2 <- as.matrix(rho2)
rotate3 <- raster(rho2[nrow(rho2):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))

am6_ns_sig <- rotate3
write.csv(rho2, file = "am6_ns_sig.csv")

###################################################################################################
# AM7
q <- 16
lonmin1 <- -94
lonmax1 <- -56

a <- readMat("SOS1_AM_R3_20N70N_7.mat")
lon <- seq(lonmin1,lonmax1,0.05)

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

a2 <- a$sos3[which(year <=2012 & year >=1981),,]
dim(a2)<- c(length(1981:2012),length(lat)*length(lon))
a2df <- as.data.frame(a2)
#test scale a2df
#a2df1 <- scale(a2df)
#a2df1df <- as.data.frame(a2df1)
#a2df <- a2df1df

m1 <- a2df[year %in% as.character(allS),]
m2 <- m1
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am7_diff_south <- rotate3
# 
m1 <- a2df[year %in% as.character(allN),]
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am7_diff_north <- rotate3
# 
diff <- colMeans(m1)-colMeans(m2)
diff <- as.matrix(diff)
diff[diff>12] <- 12
diff[-12>diff] <- -12
dim(diff) <- c(length(lat),length(lon))
diff <- as.data.frame(diff)
colnames(diff) <- lon
rownames(diff) <- lat
diff <- as.matrix(diff)
rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
am7_diff_all <- rotate3

plot(am7_diff_all)

# rho2 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))
# 
# for (i in 1:length(m1)){
#   y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T)
#   if(y < 0.1 & is.na(y) == "FALSE"){
#     rho2[i] <- mean(m1[,i])-mean(m2[,i])
#   }
# }
# #maybe the median would be better
# 
# rho1 <- as.matrix(rho2)
# dim(rho1) <- c(length(lat),length(lon))
# rho1 <- as.data.frame(rho1)
# rho1[rho1>20] <- 20
# rho1[-20>rho1] <- -20
# colnames(rho1) <- lon
# rownames(rho1) <- lat
# rho1 <- as.matrix(rho1)
# rotate3 <- raster(rho1[nrow(rho1):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# 
# am7_ns_sig <- rotate3

rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))#NA[length(mon1)]
ptm <- proc.time()   
for (i in 1:length(m1)){
  y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T) 
  z <- my.t.test.statistic(m1[,i],m2[,i],na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- z
  }
}
proc.time() - ptm

#ptm <- proc.time()   
#write.csv(rho1, file = "r7.csv")
#proc.time() - ptm

rho2 <- as.matrix(rho1)
dim(rho2) <- c(length(lat),length(lon))
rho2 <- as.data.frame(rho2)
rho2[rho2>0] <- 1
rho2[0>rho2] <- -1
colnames(rho2) <- lon
rownames(rho2) <- lat
rho2 <- as.matrix(rho2)
rotate3 <- raster(rho2[nrow(rho2):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))

am7_ns_sig <- rotate3
write.csv(rho2, file = "am7_ns_sig.csv")
#plot(am7_ns_sig)
#freq(am7_ns_sig)

# 
# ###################################################################################################
# AM8
q <- 17
lonmin1 <- -10
lonmax1 <- 8

a <- readMat("SOS1_AM_R3_20N70N_8.mat")
lon <- seq(lonmin1,lonmax1,0.05)

indices_names[q]
all <- as.data.frame(cbind(jet$YEAR,jet[,q])) #
allS <- all[order(all$V2),]$V1[1:10]
allN <- all[order(all$V2),]$V1[(length(all$V1)-10+1):length(all$V1)]

a2 <- a$sos3[which(year <=2012 & year >=1981),,]
dim(a2)<- c(length(1981:2012),length(lat)*length(lon))
a2df <- as.data.frame(a2)

m1 <- a2df[year %in% as.character(allS),]
m2 <- m1
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am8_diff_south <- rotate3
# 
m1 <- a2df[year %in% as.character(allN),]
# diff <- colMeans(m1)-colMeans(a2df)
# diff <- as.matrix(diff)
# diff[diff>12] <- 12
# diff[-12>diff] <- -12
# dim(diff) <- c(length(lat),length(lon))
# diff <- as.data.frame(diff)
# colnames(diff) <- lon
# rownames(diff) <- lat
# diff <- as.matrix(diff)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# am8_diff_north <- rotate3
# 
diff <- colMeans(m1)-colMeans(m2)
diff <- as.matrix(diff)
diff[diff>12] <- 12
diff[-12>diff] <- -12
dim(diff) <- c(length(lat),length(lon))
diff <- as.data.frame(diff)
colnames(diff) <- lon
rownames(diff) <- lat
diff <- as.matrix(diff)
rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
am8_diff_all <- rotate3


rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))#NA[length(mon1)]

for (i in 1:length(m1)){
  y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T) 
  z <- my.t.test.statistic(m1[,i],m2[,i],na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- z
  }
}

#write.csv(rho1, file = "r8.csv")

rho2 <- as.matrix(rho1)
dim(rho2) <- c(length(lat),length(lon))
rho2 <- as.data.frame(rho2)
rho2[rho2>0] <- 1
rho2[0>rho2] <- -1
colnames(rho2) <- lon
rownames(rho2) <- lat
rho2 <- as.matrix(rho2)
rotate3 <- raster(rho2[nrow(rho2):1,]) #need to flip rotate 3 on yaxis
extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))

am8_ns_sig <- rotate3
write.csv(rho2, file = "am8_ns_sig.csv")

# ## Ttest between N and S shifts
# rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))
# 
# for (i in 1:length(m1)){
#   y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T)
#   if(y < 0.1 & is.na(y) == "FALSE"){
#     rho1[i] <- mean(m1[,i])-mean(m2[,i])
#   }
# }
# #maybe the median would be better
# 
# rho1 <- as.matrix(rho1)
# dim(rho1) <- c(length(lat),length(lon))
# rho1 <- as.data.frame(rho1)
# rho1[rho1>12] <- 12
# rho1[-12>rho1] <- -12
# colnames(rho1) <- lon
# rownames(rho1) <- lat
# rho1 <- as.matrix(rho1)
# rotate3 <- raster(diff[nrow(diff):1,]) #need to flip rotate 3 on yaxis
# extent(rotate3) <- extent(c(min(lon),max(lon),min(lat),max(lat)))
# 
# am8_ns_sig <- rotate3
# 
# ##########################################################################################################################
# 
# ## Plotting difference between North and South
plot(am1_diff_all)
map("world", col = "black", add = T)
plot(am2_diff_all)
map("world", col = "black", add = T)
plot(am3_diff_all)
map("world", col = "black", add = T)
plot(am4_diff_all)
map("world", col = "black", add = T)
plot(am4b_diff_all)
map("world", col = "black", add = T)
plot(am8_diff_all)
map("world", col = "black", add = T)



# #t1 <- merge(am1,am2,am3,tolerance = 0.5)
t <- merge(am1_diff_all,am2_diff_all,am3_diff_all, tolerance = 0.5) 
plot(t)

t <- merge(am5_diff_all,am6_diff_all,am7_diff_all, tolerance = 0.5) 
# 
lon <- seq(-180,179.95,0.05)

rwb <- (brewer.pal(n = 11, name = 'RdBu'))
par(mai=c(1,1,1,1))
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F)
map.axes()


plot(am6_diff_all,col = rwb, add = T, legend = F)#,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), add = T, legend = T)

map("world", xlim=c(min(lon),max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "light gray")#, wrap=c(-180,180),add = F)
plot(t,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)#, add = T, legend = T)
plot(am4_diff_all,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)

plot(t,col = rwb, add = F, xlim=c(min(lon),max(lon)),ylim=c(min(lat),max(lat)),legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)#, add = T, legend = T)
plot(am4_diff_all,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)
map.axes()


abline(v = -150)
abline(v = -120)
abline(v = -94)
abline(v = -56)
# 
# 
# par(mai=c(0,0,0,1))  
# map("world", xlim=c(-180,180),ylim=c(20,70), fill = T, col = "black")
# map("world", xlim=c(-180,180),ylim=c(20,70), fill = T, col = "black")
# map.axes()
# 
# plot(t,col = rwb,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), add = F, legend = T)
# plot(am6_ns_sig,col = rwb,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), add = F,legend = T)
# plot(x,col = rwb,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), add = F,legend = T)
# 
# 
# for (i in 1:32){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg1)) #
#   lines(24:74,rep(all$V2[i],length(24:74)), col = "grey")
# }
# 
# for (i in 1:32){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg2)) #
#   lines(74:116,rep(all$V2[i],length(116:74)), col = "grey")
# }
# 
# for (i in 1:32){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg3)) #
#   lines(116:152,rep(all$V2[i],length(116:152)), col = "grey")
# }
# 
# for (i in 1:32){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg4)) #
#   lines(152:180,rep(all$V2[i],length(152:180)), col = "grey")
# }
# 
# for (i in 1:32){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg4)) #
#   lines(-180:-150,rep(all$V2[i],length(-180:-150)), col = "grey")
# }
# 
# for (i in 1:32){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg5)) #
#   lines(-150:-120,rep(all$V2[i],length(-150:-120)), col = "grey")
# }
# 
# for (i in 1:32){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg6)) #
#   lines(-120:-94,rep(all$V2[i],length(-120:-94)), col = "grey")
# }
# 
# for (i in 1:32){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg7)) #
#   lines(-94:-56,rep(all$V2[i],length(-94:-56)), col = "grey")
# }
# 
# for (i in 1:32){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg8)) #
#   lines(-10:8,rep(all$V2[i],length(-10:8)), col = "grey")
# }
# 
# abline(v = 24)
# abline(v = 74)
# abline(v = 116)
# abline(v = 152)
# abline(v = -150)
# abline(v = -120)
# abline(v = -94)
# abline(v = -56)
# abline(v = -10) 
# abline(v = 8) 
# 
# ###################################################################
# ## South 
# 
# t <- merge(am5_diff_south,am6_diff_south,am7_diff_south, tolerance = 0.5) 
# 
# rwb <- (brewer.pal(n = 11, name = 'RdBu'))
# 
# par(mai=c(0,0,0,1))  
# map("world", xlim=c(-180,180),ylim=c(20,70), fill = F, col = "black")
# map("world", xlim=c(-180,180),ylim=c(20,70), fill = F, col = "black", add = T)
# map.axes()
# 
# plot(t,col = rwb,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), add = F, legend = T)
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg1)) #
#   S <- all[order(all$V2),]$V2[1:10]
#   lines(24:74,rep(S[i],length(24:74)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg2)) #
#   S <- all[order(all$V2),]$V2[1:10]
#   lines(74:116,rep(S[i],length(116:74)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg3)) #
#   S <- all[order(all$V2),]$V2[1:10]
#   lines(116:152,rep(S[i],length(116:152)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg4)) #
#   S <- all[order(all$V2),]$V2[1:10]
#   lines(152:180,rep(S[i],length(152:180)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg4)) #
#   S <- all[order(all$V2),]$V2[1:10]
#   lines(-180:-150,rep(S[i],length(-180:-150)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg5)) #
#   S <- all[order(all$V2),]$V2[1:10]
#   lines(-150:-120,rep(S[i],length(-150:-120)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg6)) #
#   S <- all[order(all$V2),]$V2[1:10]
#   lines(-120:-94,rep(S[i],length(-120:-94)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg7)) #
#   S <- all[order(all$V2),]$V2[1:10]
#   lines(-94:-56,rep(S[i],length(-94:-56)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg8)) #
#   S <- all[order(all$V2),]$V2[1:10]
#   lines(-10:8,rep(S[i],length(-10:8)), col = "grey")
# }
# 
# 
# abline(v = 24)
# abline(v = 74)
# abline(v = 116)
# abline(v = 152)
# abline(v = -150)
# abline(v = -120)
# abline(v = -94)
# abline(v = -56)
# abline(v = -10) 
# abline(v = 8) 
# 
# map.axes()
# 
# ############################################################################################################### North
# 
# t <- merge(am5_diff_north,am6_diff_north,am7_diff_north, tolerance = 0.5) 
# 
# rwb <- (brewer.pal(n = 11, name = 'RdBu'))
# 
# par(mai=c(0,0,0,1))  
# map("world", xlim=c(-180,180),ylim=c(20,70), fill = F, col = "black")
# map("world", xlim=c(-180,180),ylim=c(20,70), fill = F, col = "black", add = T)
# map.axes()
# 
# plot(t,col = rwb,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), add = F, legend = T)
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg1)) #
#   S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
#   lines(24:74,rep(S[i],length(24:74)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg2)) #
#   S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
#   lines(74:116,rep(S[i],length(116:74)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg3)) #
#   S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
#   lines(116:152,rep(S[i],length(116:152)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg4)) #
#   S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
#   lines(152:180,rep(S[i],length(152:180)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg4)) #
#   S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
#   lines(-180:-150,rep(S[i],length(-180:-150)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg5)) #
#   S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
#   lines(-150:-120,rep(S[i],length(-150:-120)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg6)) #
#   S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
#   lines(-120:-94,rep(S[i],length(-120:-94)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg7)) #
#   S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
#   lines(-94:-56,rep(S[i],length(-94:-56)), col = "grey")
# }
# 
# for (i in 1:10){
#   all <- as.data.frame(cbind(jet$YEAR,jet$AM_Reg8)) #
#   S <- all[order(all$V2),]$V2[(length(all$V2)-10+1):length(all$V2)]
#   lines(-10:8,rep(S[i],length(-10:8)), col = "grey")
# }
# 
# 
# abline(v = 24)
# abline(v = 74)
# abline(v = 116)
# abline(v = 152)
# abline(v = -150)
# abline(v = -120)
# abline(v = -94)
# abline(v = -56)
# abline(v = -10) 
# abline(v = 8) 
# 
# map.axes()
# 

## SIG 
plot(am1_ns_sig)
map("world", col = "black", add = T)
plot(am2_ns_sig)
map("world", col = "black", add = T)
plot(am3_ns_sig)
map("world", col = "black", add = T)
plot(am4a_ns_sig)
map("world", col = "black", add = T)
plot(am4b_ns_sig)
map("world", col = "black", add = T)



# #t1 <- merge(am1,am2,am3,tolerance = 0.5)
t <- merge(am1_ns_sig,am2_ns_sig,am3_ns_sig, tolerance = 0.5) 
plot(t)

t <- merge(am5_diff_all,am6_diff_all,am7_diff_all, tolerance = 0.5) 
# 
lon <- seq(-180,179.95,0.05)

rwb <- (brewer.pal(n = 11, name = 'RdBu'))
par(mai=c(1,1,1,1))
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F)
map.axes()


plot(am6_diff_all,col = rwb, add = T, legend = F)#,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), add = T, legend = T)

map("world", xlim=c(min(lon),max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "light gray")#, wrap=c(-180,180),add = F)
plot(t,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)#, add = T, legend = T)
plot(am4_diff_all,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)

plot(t,col = rwb, add = F, xlim=c(min(lon),max(lon)),ylim=c(min(lat),max(lat)),legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)#, add = T, legend = T)
plot(am4_diff_all,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)
map.axes()




## experiment to try and get all rasters in 1 file

library(gdalUtils)
library(rgdal)
all_my_rasts <- c('am1_diff_all', 'am2_diff_all', 'am3_diff_all','am4_diff_all')
e <- extent(0, 180, 20, 70)
template <- raster(e)
projection(template) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
writeRaster(template, file="amsos1.tif", format="GTiff")
mosaic_rasters(gdalfile=all_my_rasts,dst_dataset="amsos1.tif",of="GTiff")
gdalinfo("amsos1.tif")
#Reduce(function(...)merge(...,tolerance=1),ast14dmo.sd)

png("sos1am_a.png",13,8,
    units = "in",res = 600, pointsize=20, family= "helvetica")
#par(mfrow=c(1,2), tcl=-0.5, family="serif", mai=c(0,0,0,0),mar = c(0, 0, 0, 0))
map("world", xlim=c(-180,180),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = F)#, wrap=c(-180,180),add = F)
map.axes()
#t <- merge(am1_diff_all,am2_diff_all,am3_diff_all,am4_diff_all, tolerance = 1) 
#plot(t)
#map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)
t <- merge(am1_diff_all,am2_diff_all,am3_diff_all, tolerance = 1) 
plot(t, add = T, col = rwb, legend = F)
plot(am4_diff_all, add = T, col = rwb, legend = F)
dev.off()

#maybe set rasters to have all the same extent? but would that warp them?
#setEPS()
#postscript("sos1am_b.eps")

png("sos1am_b.png",13,8,
    units = "in",res = 600, pointsize=20, family= "helvetica")
#par(mfrow=c(1,2), tcl=-0.5, family="serif", mai=c(0,0,0,0),mar = c(0, 0, 0, 0))
map("world", xlim=c(-180,180),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = F)#, wrap=c(-180,180),add = F)
map.axes()
#t <- merge(am1_diff_all,am2_diff_all,am3_diff_all,am4_diff_all, tolerance = 1) 
#plot(t)
#map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)
t <- merge(am4b_diff_all,am5_diff_all,am6_diff_all, am7_diff_all, tolerance = 1) 
plot(t, add = T, col = rwb,legend = F)
plot(am8_diff_all, add = T,col = rwb,legend = F)
#map("world", xlim=c(-180,180),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)

#abline(v = -150)
#abline(v = -120)
#abline(v = -94)
#abline(v = -56)
dev.off()


#rm(list = ls(all.names = TRUE))
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
#ON_1 
q <- 26
lonmin1 <- 10
lonmax1 <- 26
a <- readMat("EOS1_R3_20N70N_ON_1.mat")
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
# on1_diff_south <- rotate3

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
# on1_diff_north <- rotate3

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
on1_diff_all <- rotate3

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

on1_ns_sig <- rotate3
write.csv(rho2, file = "on1_ns_sig.csv")
# ##############################################################
# ###################################################################################################
# # NOW AM2
q <- 27
lonmin1 <- 26
lonmax1 <- 56
a <- readMat("EOS1_R3_20N70N_ON_2.mat")
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
# on2_diff_south <- rotate3
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
# on2_diff_north <- rotate3
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
on2_diff_all <- rotate3

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

on2_ns_sig <- rotate3
write.csv(rho2, file = "on2_ns_sig.csv")
# ###################################################################################################
# AM3

q <- 28
lonmin1 <- 56
lonmax1 <- 104
a <- readMat("EOS1_R3_20N70N_ON_3.mat")
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
# on3_diff_south <- rotate3
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
# on3_diff_north <- rotate3
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
on3_diff_all <- rotate3
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

on3_ns_sig <- rotate3
write.csv(rho2, file = "on3_ns_sig.csv")
# ###################################################################################################
# AM4 #part 1
q <- 29
lonmin1 <- 104
lonmax1 <- 148

a <- readMat("EOS1_R3_20N70N_ON_4.mat")
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
# on4_diff_south <- rotate3
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
# on4_diff_north <- rotate3
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
on4_diff_all <- rotate3

plot(on4_diff_all)

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

on4_ns_sig <- rotate3
plot(on4a_ns_sig, col = rwb)
plot(on4_diff_all, col = rwb)

write.csv(rho2, file = "on4_ns_sig.csv")

# ############################################################################
# ###################################################################################################
# ON5
q <- 30
lonmin1 <- 148
lonmax1 <- 180

a <- readMat("EOS1_R3_20N70N_ON_5.mat")
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
# on4b_diff_south <- rotate3
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
# on4b_diff_north <- rotate3
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
on5_diff_all <- rotate3

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

on5_ns_sig <- rotate3
write.csv(rho2, file = "on5_ns_sig.csv")
# ###################################################################################################
# ON5b
q <- 30
lonmin1 <- -180
lonmax1 <- -146

a <- readMat("EOS1_R3_20N70N_ON_5b.mat") #
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
on5b_diff_all <- rotate3

plot(on5b_diff_all)
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

on5b_ns_sig <- rotate3
write.csv(rho2, file = "on5b_ns_sig.csv")
# 
# ###################################################################################################
# 
# ON6
q <- 31
lonmin1 <- -146
lonmax1 <- -98
a <- readMat("EOS1_R3_20N70N_ON_6.mat")
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
on6_diff_all <- rotate3

plot(on6_diff_all)
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

on6_ns_sig <- rotate3
write.csv(rho2, file = "on6_ns_sig.csv")

###################################################################################################
# ON7
q <- 32
lonmin1 <- -98
lonmax1 <- -66

a <- readMat("EOS1_R3_20N70N_ON_7.mat")
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
on7_diff_all <- rotate3

plot(on7_diff_all)

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

on7_ns_sig <- rotate3
write.csv(rho2, file = "on7_ns_sig.csv")
#plot(am7_ns_sig)
#freq(am7_ns_sig)

# 
# ###################################################################################################
# AM8
q <- 33
lonmin1 <- -66
lonmax1 <- -34

a <- readMat("EOS1_R3_20N70N_ON_8.mat")
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
on8_diff_all <- rotate3


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

on8_ns_sig <- rotate3
write.csv(rho2, file = "on8_ns_sig.csv")

# ## Ttest between N and S shifts
# rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))
# 
# for (i in 1:length(m1)){
#   y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T)
#   if(y < 0.1 & is.na(y) == "FALSE"){
#     rho1[i] <- mean(m1[,i])-mean(m2[,i])
#   }
# }
# # ###################################################################################################
# ON9
q <- 34
lonmin1 <- -34
lonmax1 <- -10

a <- readMat("EOS1_R3_20N70N_ON_8.mat")
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
on9_diff_all <- rotate3


rho1 <- as.data.frame(matrix(NA,nrow = 1,ncol = length(m1)))#NA[length(mon1)]

for (i in 1:length(m1)){
  y <- my.t.test.p.value(m1[,i],m2[,i],na.omit = T) 
  z <- my.t.test.statistic(m1[,i],m2[,i],na.omit = T)
  if(y <= 0.1 & is.na(y) == "FALSE"){
    rho1[i] <- z
  }
}

#write.csv(rho1, file = "r9.csv")

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

on9_ns_sig <- rotate3
write.csv(rho2, file = "on9_ns_sig.csv")

# ##########################################################################################################################
# 
# ## Plotting difference between North and South
plot(on1_diff_all)
map("world", col = "black", add = T)
plot(on2_diff_all)
map("world", col = "black", add = T)
plot(on3_diff_all)
map("world", col = "black", add = T)
plot(on4_diff_all)
map("world", col = "black", add = T)
plot(on5_diff_all)
map("world", col = "black", add = T)
plot(on6_diff_all)
map("world", col = "black", add = T)
plot(on7_diff_all)
map("world", col = "black", add = T)
plot(on8_diff_all)
map("world", col = "black", add = T)
plot(on9_diff_all)
map("world", col = "black", add = T)

# #t1 <- merge(on1,on2,on3,tolerance = 0.5)
t <- merge(on1_diff_all,on2_diff_all,on3_diff_all, tolerance = 0.5) 
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
plot(on4_diff_all,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)

plot(t,col = rwb, add = F, xlim=c(min(lon),max(lon)),ylim=c(min(lat),max(lat)),legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)#, add = T, legend = T)
plot(on4_diff_all,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)
map.axes()


abline(v = -150)
abline(v = -120)
abline(v = -94)
abline(v = -56)
# 
# 

## SIG 
plot(on1_ns_sig)
map("world", col = "black", add = T)
plot(on2_ns_sig)
map("world", col = "black", add = T)
plot(on3_ns_sig)
map("world", col = "black", add = T)
plot(on4a_ns_sig)
map("world", col = "black", add = T)
plot(on4b_ns_sig)
map("world", col = "black", add = T)



# #t1 <- merge(on1,on2,on3,tolerance = 0.5)
t <- merge(on1_ns_sig,on2_ns_sig,on3_ns_sig, tolerance = 0.5) 
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
plot(on4_diff_all,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)

plot(t,col = rwb, add = F, xlim=c(min(lon),max(lon)),ylim=c(min(lat),max(lat)),legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)#, add = T, legend = T)
plot(on4_diff_all,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)
map.axes()


png("eos1on_a.png",13,8,
    units = "in",res = 600, pointsize=20, family= "helvetica")
#par(mfrow=c(1,2), tcl=-0.5, family="serif", mai=c(0,0,0,0),mar = c(0, 0, 0, 0))
map("world", xlim=c(-180,180),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = F)#, wrap=c(-180,180),add = F)
map.axes()
#t <- merge(on1_diff_all,on2_diff_all,on3_diff_all,on4_diff_all, tolerance = 1) 
#plot(t)
#map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)
t <- merge(on1_diff_all,on2_diff_all,on3_diff_all,on4_diff_all, tolerance = 1) 
plot(t, add = T, col = rwb, legend = F)
plot(on5_diff_all, add = T, col = rwb, legend = F)
dev.off()

#maybe set rasters to have all the same extent? but would that warp them?
#setEPS()
#postscript("sos1am_b.eps")

png("eos1on_b.png",13,8,
    units = "in",res = 600, pointsize=20, family= "helvetica")
#par(mfrow=c(1,2), tcl=-0.5, family="serif", mai=c(0,0,0,0),mar = c(0, 0, 0, 0))
map("world", xlim=c(-180,180),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = F)#, wrap=c(-180,180),add = F)
map.axes()
#t <- merge(on1_diff_all,on2_diff_all,on3_diff_all,on4_diff_all, tolerance = 1) 
#plot(t)
#map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)
t <- merge(on5b_diff_all,on6_diff_all, on7_diff_all,on8_diff_all, tolerance = 1) 
plot(t, add = T, col = rwb,legend = F)
plot(on9_diff_all, add = T,col = rwb,legend = F)
#map("world", xlim=c(-180,180),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)

#abline(v = -150)
#abline(v = -120)
#abline(v = -94)
#abline(v = -56)
dev.off()


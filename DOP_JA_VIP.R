
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
#JA_1 
q <- 19
lonmin1 <- 20
lonmax1 <- 53.95 #not sure why this is happening... difference between this and sos? #check matlab?
a <- readMat("DOP1_R3_20N70N_JA_1.mat")
lon <- seq(lonmin1,lonmax1,0.05) #179.95 lon

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
# ja1_diff_south <- rotate3

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
ja1_diff_all <- rotate3

plot(ja1_diff_all)
map("world", xlim=c(min(lon),max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)



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

ja1_ns_sig <- rotate3
write.csv(rho2, file = "ja1_ns_sig.csv")
# ##############################################################
# ###################################################################################################
# # NOW JA2
q <- 20
lonmin1 <- 54
lonmax1 <- 85.95
a <- readMat("DOP1_R3_20N70N_JA_2.mat")
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
ja2_diff_all <- rotate3

plot(ja2_diff_all)
map("world", xlim=c(min(lon),max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)

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

ja2_ns_sig <- rotate3
write.csv(rho2, file = "ja2_ns_sig.csv")
# ###################################################################################################
# JA3

q <- 21
lonmin1 <- 86
lonmax1 <- 180
a <- readMat("DOP1_R3_20N70N_JA_3.mat")
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
ja3_diff_all <- rotate3
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

ja3_ns_sig <- rotate3
write.csv(rho2, file = "ja3_ns_sig.csv")
# ###################################################################################################
# JA3 B
q <- 21
lonmin1 <- -180
lonmax1 <- -160

a <- readMat("DOP1_R3_20N70N_JA_3b.mat")
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
# ja3b_diff_north <- rotate3
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
ja3b_diff_all <- rotate3

plot(ja3b_diff_all)

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

ja3b_ns_sig <- rotate3
plot(ja3b_ns_sig, col = rwb)
plot(ja3b_diff_all, col = rwb)

write.csv(rho2, file = "ja3b_ns_sig.csv")

# ############################################################################
# ###################################################################################################
# JA4 
q <- 22
lonmin1 <- -160
lonmax1 <- -104

a <- readMat("DOP1_R3_20N70N_JA_4.mat")
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
ja4_diff_all <- rotate3

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

ja4_ns_sig <- rotate3
write.csv(rho2, file = "ja4_ns_sig.csv")
# ###################################################################################################
# JA5
q <- 23
lonmin1 <- -104
lonmax1 <- -58

a <- readMat("DOP1_R3_20N70N_JA_5.mat") #
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
ja5_diff_all <- rotate3

plot(ja5_diff_all)
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

ja5_ns_sig <- rotate3
write.csv(rho2, file = "ja5_ns_sig.csv")
# 
# ###################################################################################################
# 
# JA6
q <- 24
lonmin1 <- -58
lonmax1 <- -16
a <- readMat("DOP1_R3_20N70N_JA_6.mat")
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
ja6_diff_all <- rotate3

plot(ja6_diff_all)
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

ja6_ns_sig <- rotate3
write.csv(rho2, file = "ja6_ns_sig.csv")

###################################################################################################
# JA7
q <- 25
lonmin1 <- -16
lonmax1 <- 8

a <- readMat("DOP1_R3_20N70N_JA_7.mat")
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
ja7_diff_all <- rotate3

plot(ja7_diff_all)

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

ja7_ns_sig <- rotate3
write.csv(rho2, file = "ja7_ns_sig.csv")
#plot(ja7_ns_sig)
#freq(ja7_ns_sig)

############################################################################################### 

## SIG 
# plot(am1_ns_sig)
# map("world", col = "black", add = T)
# plot(am2_ns_sig)
# map("world", col = "black", add = T)
# plot(am3_ns_sig)
# map("world", col = "black", add = T)
# plot(ja3b_ns_sig)
# map("world", col = "black", add = T)
# plot(am4b_ns_sig)
# map("world", col = "black", add = T)



# #t1 <- merge(am1,am2,am3,tolerance = 0.5)
#t <- merge(am1_ns_sig,am2_ns_sig,am3_ns_sig, tolerance = 0.5) 
#plot(t)

t <- merge(ja5_diff_all,ja6_diff_all,ja7_diff_all, tolerance = 0.5) 
# 
lon <- seq(-180,179.95,0.05)

rwb <- (brewer.pal(n = 11, name = 'RdBu'))
par(mai=c(1,1,1,1))
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = TRUE, col = "light gray", border="light gray")#, wrap=c(-180,180),add = F)
map.axes()


plot(ja6_diff_all,col = rwb, add = T, legend = F)#,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), add = T, legend = T)

map("world", xlim=c(min(lon),max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "light gray")#, wrap=c(-180,180),add = F)
plot(t,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)#, add = T, legend = T)
plot(ja4_diff_all,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)

plot(t,col = rwb, add = F, xlim=c(min(lon),max(lon)),ylim=c(min(lat),max(lat)),legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)#, add = T, legend = T)
plot(ja4_diff_all,col = rwb, add = T, legend = F,breaks = c(-12,-10,-8,-6,-4,-2,2,4,6,8,10,12), axes = F)
map("world", xlim=c(0,max(lon)),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = T)#, wrap=c(-180,180),add = F)
map.axes()





png("dop1ja_a.png",13,8,
    units = "in",res = 600, pointsize=20, family= "helvetica")
map("world", xlim=c(-180,180),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = F)#, wrap=c(-180,180),add = F)
map.axes()
t <- merge(ja1_diff_all,ja2_diff_all,ja3_diff_all, tolerance = 1) 
plot(t, add = T, col = rwb, legend = F)
dev.off()

#maybe set rasters to have all the same extent? but would that warp them?
#setEPS()
#postscript("sos1am_b.eps")

png("dop1ja_b.png",13,8,
    units = "in",res = 600, pointsize=20, family= "helvetica")
map("world", xlim=c(-180,180),ylim=c(min(lat),max(lat)), fill = F, col = "black", add = F)#, wrap=c(-180,180),add = F)
map.axes()
t <- merge(ja4_diff_all,ja5_diff_all,ja6_diff_all, ja7_diff_all, tolerance = 1) 
plot(t, add = T, col = rwb,legend = F)
plot(ja3b_diff_all, add = T,col = rwb,legend = F)
dev.off()


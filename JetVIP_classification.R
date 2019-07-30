# Land Cover, Koppen-Greiger climate classification, and Elevation test

# So Far for AM6 only__ modifying for all___ 1/24/19 3:30pm

# I've subsetted the MODIS LC into AM, JA, and JF regions in matlab, saved as .mat files, and now read them into R

#Import the results of the significance file -1 and 1

#MAY BE USEFUL TO HAVE A EARLY LATE %TOTAL BY REGION for LC and KG
#what % of significant pixels are early or late?

f <- system.file("external/test.grd", package="raster")
r <- raster(f)
x <- rasterToContour(r)
class(x)
plot(r)
plot(x, add=TRUE)


am6 <- read.csv("on7_ns_sig.csv")
am6[,1] <- NULL

am6 <- read.csv("am5_ns_sig.csv")
am6[,1] <- NULL

am6 <- read.csv("am6_ns_sig.csv")
am6[,1] <- NULL

am6 <- as.data.frame(as.matrix(am1_ns_sig))
am6 <- as.data.frame(as.matrix(am2_ns_sig))
am6 <- as.data.frame(as.matrix(am3_ns_sig))
am6 <- as.data.frame(as.matrix(am4a_ns_sig))
am6 <- as.data.frame(as.matrix(am4b_ns_sig))
am6 <- as.data.frame(as.matrix(am7_ns_sig))
am6 <- as.data.frame(as.matrix(am8_ns_sig))



# subsetting am6 by pixels that have an average start of spring in AM
# library(R.matlab)
# a <- readMat("SOS1_R3_20N70N_6.mat")
# lonmin1 <- -120
# lonmax1 <- -94
# lon <- seq(lonmin1,lonmax1,0.05)
# year <- 1981:2014
# lat <- seq(20,70,0.05)
# 
# a2 <- a$sos3[which(year <=2012 & year >=1981),,]
# dim(a2)<- c(length(1981:2012),length(lat)*length(lon))
# a2df <- as.data.frame(a2)
# colMeans(a2)
# a2_colmeans <- colMeans(a2,na.rm = T)
# 
# a2df1 <- a2df
# 
# am <- which(a2_colmeans<151 & a2_colmeans>91)
# dim(a2)
# 
# am2NA <- matrix(NA,nrow = 1,ncol = 521521)
# dim(am2NA)
# am2NA[,which(a2_colmeans<151 & a2_colmeans>91)] <- 1
# dim(am2NA) <- c(length(lat),length(lon)) #1 and NA
# 
# b <- am6*am2NA
# freq(raster(as.matrix(b)))
# freq(raster(as.matrix(am6)))
# 
# am6 <- b
###############################################################################
library(R.matlab)
a <- readMat("MODISLC_AM_1.mat")# proj=GCTP_GEO
a <- readMat("MODISLC_AM_2.mat")# proj=GCTP_GEO
a <- readMat("MODISLC_AM_3.mat")# proj=GCTP_GEO
a <- readMat("MODISLC_AM_4.mat")# proj=GCTP_GEO
a <- readMat("MODISLC_AM_4b.mat")# proj=GCTP_GEO
a <- readMat("MODISLC_AM_5.mat")# proj=GCTP_GEO
a <- readMat("MODISLC_AM_6.mat")# proj=GCTP_GEO
a <- readMat("MODISLC_AM_7.mat")# proj=GCTP_GEO
a <- readMat("MODISLC_AM_8.mat")# proj=GCTP_GEO

a <- a$land2

key <- as.data.frame(cbind(0:16, c("water","evergreen needleleaf forest","evergreen broadleaf forest","deciduous needleleaf forest","deciduous broadleaf forest","mixed forests","closed shrubland","open shrublands","woody savannas","savannas","grasslands","permanent wetlands","croplands","urban and built-up","cropland/natural vegetation mosaic","snow and ice","barren or sparsely vegetated")))
colnames(key) <- c('SYM','LC')

den <- as.data.frame(freq(raster(as.matrix(a))))

key[key$SYM %in% den$value,3] <- den$count

dim(a)
b <- am6*a
#freq(raster(as.matrix(am6)))
#freq(raster(as.matrix(a)))

freq(raster(as.matrix(b)))

# match the value column of the numerator with the value of the denominator, divide count of that num by count of denom

num_pos <- as.data.frame(freq(raster(as.matrix(b)))[which(freq(raster(as.matrix(b)))[,1]>0),])
lc_pos <- as.data.frame(num_pos$count/den[den$value %in% num_pos$value,]$count*100)
rownames(lc_pos) <- den[den$value %in% num_pos$value,]$value
lc_pos[,1] <- round(lc_pos[,1],2)

num_neg <- as.data.frame(freq(raster(as.matrix(b)))[which(freq(raster(as.matrix(b)))[,1]<0),])
num_neg[,1] <- abs(num_neg[,1])
num_neg <- num_neg[order(num_neg$value),]
#good till here
lc_neg <- as.data.frame(num_neg$count/den[den$value %in% num_neg$value,]$count*100)
rownames(lc_neg) <- den[den$value %in% num_neg$value,]$value
lc_neg[,1] <- round(lc_neg[,1],2)

#merge these stats into one table
 
#how to fill in the data frame with missing sequential values? will that help me then match/merge them and divide by the count?

key[key$SYM %in% num_neg$value,4] <- lc_neg[,1]
key[key$SYM %in% num_pos$value,5] <- lc_pos[,1]

colnames(key) <- c('SYM','LC_TYPE','AM6_TOTAL_PIXELS','AM6_DIFF_%EARLY','AM6_DIFF_%LATE')
colnames(key) <- c('SYM','LC_TYPE','AM1_TOTAL_PIXELS','AM1_DIFF_%EARLY','AM1_DIFF_%LATE')
colnames(key) <- c('SYM','LC_TYPE','AM2_TOTAL_PIXELS','AM2_DIFF_%EARLY','AM2_DIFF_%LATE')
colnames(key) <- c('SYM','LC_TYPE','AM3_TOTAL_PIXELS','AM3_DIFF_%EARLY','AM3_DIFF_%LATE')
colnames(key) <- c('SYM','LC_TYPE','AM4a_TOTAL_PIXELS','AM4a_DIFF_%EARLY','AM4a_DIFF_%LATE')
colnames(key) <- c('SYM','LC_TYPE','AM4b_TOTAL_PIXELS','AM4b_DIFF_%EARLY','AM4b_DIFF_%LATE')
colnames(key) <- c('SYM','LC_TYPE','AM5_TOTAL_PIXELS','AM5_DIFF_%EARLY','AM5_DIFF_%LATE')
colnames(key) <- c('SYM','LC_TYPE','AM7_TOTAL_PIXELS','AM7_DIFF_%EARLY','AM7_DIFF_%LATE')
colnames(key) <- c('SYM','LC_TYPE','AM8_TOTAL_PIXELS','AM8_DIFF_%EARLY','AM8_DIFF_%LATE')


#work on adding plus one column each time? may not be worth the coding hours
key[,3:5]

#am8_diff_LC <- key
#am7_diff_LC <- key
#am6_diff_LC <- key
#am5_diff_LC <- key
#am4b_diff_LC <- key
#am4a_diff_LC <- key
#am3_diff_LC <- key
#am2_diff_LC <- key
#am1_diff_LC <- key

#LC <- cbind(am1_diff_LC, am2_diff_LC[,3:5],am3_diff_LC[,3:5],am4a_diff_LC[,3:5], am4b_diff_LC[,3:5],am5_diff_LC[,3:5],am6_diff_LC[,3:5], am7_diff_LC[,3:5],am8_diff_LC[,3:5])


# # elevation test
# setwd("/Volumes/AOP-NEON1.4/VIP/")
# setwd("~/Downloads/")
# 
# e <- read.csv("e10g")
# 
# install.packages("elevatr")
# library(elevatr)
# 
# data(lake)
# plot(get_elev_raster(lake, z = 2))
# examp_df <- data.frame(x = runif(10, min = -73, max = -71), y = runif(10, min = 41,max = 45))
# prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# plot(get_elev_raster(examp_df,prj = prj_dd,z =5))
# 
# install.packages("dismo")
# library(dismo)
# 
# setwd("/Volumes/AOP-NEON1.4/VIP/JetVIP2/")
# 
# library(raster)
# getData('ISO3')
# 
# #future problem of europe having so many countries for getting data:
# # library(raster)
# # misc = list()
# # misc$countries = c("ZAF", "LSO", "SWZ", "ZWE", "MOZ", "NAM", "BWA")
# # ctry_shps = do.call("bind", lapply(misc$countries, 
# #                                    function(x) getData('GADM', country=x, level=0)))
# 
# library(dplyr)
# library(purrr)
# 
# 
# elevation <- getData("alt", country = "CAN")
# 
# #x <- terrain(elevation, opt = c("slope", "aspect"), unit = "degrees")
# plot(elevation)
# 
# elevation1 <- getData("alt", country = "USA")
# elevation_a <- elevation1$`/Users/amyhudson/Documents/JetVIP/USA1_msk_alt.grd`
# elevation_b <- elevation1$`/Users/amyhudson/Documents/JetVIP/USA2_msk_alt.grd`
# elevation_c <- elevation1$`/Users/amyhudson/Documents/JetVIP/USA3_msk_alt.grd`
# elevation_d <- elevation1$`/Users/amyhudson/Documents/JetVIP/USA4_msk_alt.grd`
# 
# elevation_a <- elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP/USA1_msk_alt.grd`
# elevation_b <- elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP/USA2_msk_alt.grd`
# elevation_c <- elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP/USA3_msk_alt.grd`
# elevation_d <- elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP/USA4_msk_alt.grd`
# #plot(elevation_a, add = T) #but different color scheme, so need to change for plotting on same map
# #plot(elevation_b, add = T)
# 
# elevation2 <- getData("alt", country = "MEX")
# 
# a <- merge(elevation,elevation_a,elevation_b,elevation_c,elevation_d,elevation2)
# 
# a@extent
# 
# am6 <- as.matrix(am6)
# am6rotate <- raster(am6[nrow(am6):1,]) #need to flip rotate 3 on yaxis
# e <- extent(-120,-94,20,70)
# e <- extent(74,116,20,70) 
# extent(am6rotate) <- e
# #crop extent
# lonmin1 <- -120
# lonmax1 <- -94
# newext <- c(lonmin1, lonmax1, 20, 70)
# 
# am6_elevation <- crop(a,newext)
# 
# am6_elevation@ncols
# am6_elevation@nrows
# 
# crs(am6_elevation)
# plot(am6_elevation)
# 
# am6_elev_resampled <- resample(am6_elevation,am6rotate) #bilinear interpolation 
# am6_elev_resampled@extent
# 
# b <- as.matrix(am6_elev_resampled) * am6
# 
# b <- as.matrix(am6_elev_resampled) * abs(am6) #run this for plotting elevation by doy
# summary(b)
# as.matrix(b)
# 
# c <- as.numeric(unlist(b))
# c <- c[complete.cases(c)]
# boxplot(c, main = "N v S Significantly early/ late multiplied by elevation")
# 
# diff <- colMeans(m1)-colMeans(a2df)
# d <- as.numeric(unlist(diff))
# d[d == "NaN"] = "NA"  
# d <- as.numeric(d)
# cd <- data.frame(cbind(c,d))
# cd_narm <- cd[complete.cases(cd),]
# plot(cd_narm$d,cd_narm$c, xlab = "DOY",ylab = "Elevation", main = "AM Mask")
# legend("topright",
#        c("colMeans(N)-colMeans(ALL)","colMeans(S)-colMeans(ALL)"),
#        pch=c(1,1), # gives the legend appropriate symbols (lines)
#        col=c("black","red"),
#        cex = 0.75) # gives the legend lines the correct color and width
# 
# diff <- colMeans(m2)-colMeans(a2df)
# d <- as.numeric(unlist(diff))
# d[d == "NaN"] = "NA"  
# d <- as.numeric(d)
# cd <- data.frame(cbind(c,d))
# cd_narm <- cd[complete.cases(cd),]
# points(cd_narm$d,cd_narm$c, col = "red")
# 
# #length(seq(20,70,0.05)) #1001
# #length(seq(-120,-94,0.05)) #521
# 
# spTransform(am6_elevation)
# ##############################################################################################################
## Koppen Geiger
#Beck, H.E., N.E. Zimmermann, T.R. McVicar, N. Vergopolan, A. Berg, E.F. Wood:Present and future Köppen-Geiger climate classification maps at 1-km resolution,Nature Scientific Data, 2018.
#install.packages("rgdal")
library(rgdal)
library(raster)
#t083<- raster('Beck_KG_V1/Beck_KG_V1_present_0p083.tif')
#plot(t083) This only has 30 
#cc <- t083

cc <- raster('KOPPEN_GEO_05.tif') 
# I think Bill gave me a dated version of KG because of the 31 32 layers
# Peel MC, Finlayson BL & McMahon TA (2007), Updated world map of the Köppen-Geiger climate classification, Hydrol. Earth Syst. Sci., 11, 1633-1644.
#https://people.eng.unimelb.edu.au/mpeel/koppen.html
extent(cc)
e <- extent(-120,-94,20,70) #AM6
e <- extent(24,74,20,70) #AM1
e <- extent(74,116,20,70) #AM2
e <- extent(116,152,20,70) #AM3
e <- extent(152,180,20,70) #AM4a
e <- extent(-180,-150,20,70) #AM4b
e <- extent(-150,-120,20,70) #AM5
e <- extent(-120,-94,20,70) #AM6
e <- extent(-94,-56,20,70) #AM7
e <- extent(-10,8,20,70) #AM8

cccropped <- crop(cc,e) #ncols 520 nrows 1000 
cccropped@data@values[cccropped@data@values==255] <- 0

#freq(raster(as.matrix(cc)))
#freq(raster(as.matrix(cccropped)))

key <- as.data.frame(cbind(0:32, c("NA","Af","Am","Aw/As","BWh","BWk","BSh","BSk","Csa","Csb","Csc","Cwa","Cwb","Cwc","Cfa","Cfb","Cfc","Dsa","Dsb","Dsc","Dsd","Dwa","Dwb","Dwc","Dwd","Dfa","Dfb","Dfc","Dfd","ET","EF","ET>1500m","EF>1500m")))
colnames(key) <- c('SYM','KG')

den <- as.data.frame(freq(raster(as.matrix(cccropped))))

key[key$SYM %in% den$value,3] <- den$count

dim(cccropped)
b <- am6*cccropped
#freq(raster(as.matrix(b)))

# match the value column of the numerator with the value of the denominator, divide count of that num by count of denom

num_pos <- as.data.frame(freq(raster(as.matrix(b)))[which(freq(raster(as.matrix(b)))[,1]>0),])
lc_pos <- as.data.frame(num_pos$count/den[den$value %in% num_pos$value,]$count*100)
rownames(lc_pos) <- den[den$value %in% num_pos$value,]$value
lc_pos[,1] <- round(lc_pos[,1],2)

num_neg <- as.data.frame(freq(raster(as.matrix(b)))[which(freq(raster(as.matrix(b)))[,1]<0),])
num_neg[,1] <- abs(num_neg[,1])
num_neg <- num_neg[order(num_neg$value),]

lc_neg <- as.data.frame(num_neg$count/den[den$value %in% num_neg$value,]$count*100)
rownames(lc_neg) <- den[den$value %in% num_neg$value,]$value
lc_neg[,1] <- round(lc_neg[,1],2)

#merge these stats into one table

#how to fill in the data frame with missing sequential values? will that help me then match/merge them and divide by the count?

key[key$SYM %in% num_neg$value,4] <- lc_neg[,1]
key[key$SYM %in% num_pos$value,5] <- lc_pos[,1]

colnames(key) <- c('SYM','KG_TYPE','AM6_TOTAL_PIXELS','AM6_DIFF_%EARLY','AM6_DIFF_%LATE')
colnames(key) <- c('SYM','KG_TYPE','AM1_TOTAL_PIXELS','AM1_DIFF_%EARLY','AM1_DIFF_%LATE')
colnames(key) <- c('SYM','KG_TYPE','AM2_TOTAL_PIXELS','AM2_DIFF_%EARLY','AM2_DIFF_%LATE')
#work on adding plus one column each time? may not be worth the coding hours
colnames(key) <- c('SYM','KG_TYPE','AM3_TOTAL_PIXELS','AM3_DIFF_%EARLY','AM3_DIFF_%LATE')
colnames(key) <- c('SYM','KG_TYPE','AM4a_TOTAL_PIXELS','AM4a_DIFF_%EARLY','AM4a_DIFF_%LATE')
colnames(key) <- c('SYM','KG_TYPE','AM4b_TOTAL_PIXELS','AM4b_DIFF_%EARLY','AM4b_DIFF_%LATE')
colnames(key) <- c('SYM','KG_TYPE','AM5_TOTAL_PIXELS','AM5_DIFF_%EARLY','AM5_DIFF_%LATE')
colnames(key) <- c('SYM','KG_TYPE','AM7_TOTAL_PIXELS','AM7_DIFF_%EARLY','AM7_DIFF_%LATE')
colnames(key) <- c('SYM','KG_TYPE','AM8_TOTAL_PIXELS','AM8_DIFF_%EARLY','AM8_DIFF_%LATE')


#work on adding plus one column each time? may not be worth the coding hours
#am8_diff_KG <- key
#am7_diff_KG <- key
#am6_diff_KG <- key
#am5_diff_KG <- key
#am4b_diff_KG <- key
#am4a_diff_KG <- key
#am3_diff_KG <- key
#am2_diff_KG <- key
#am1_diff_KG <- key

KG <- cbind(am1_diff_KG, am2_diff_KG[,3:5],am3_diff_KG[,3:5],am4a_diff_KG[,3:5], am4b_diff_KG[,3:5],am5_diff_KG[,3:5],am6_diff_KG[,3:5], am7_diff_KG[,3:5],am8_diff_KG[,3:5])

# Land Cover, Koppen-Greiger climate classification, and Elevation test

# I've subsetted the MODIS LC into AM, JA, and JF regions in matlab, saved as .mat files, and now read them into R

#Import the results of the significance file -1 and 1

am6 <- read.csv("am6_ns_sig.csv")
am6[,1] <- NULL

library(R.matlab)
a <- readMat("MODISLC_AM_6.mat")# proj=GCTP_GEO
a <- a$land2

key <- as.data.frame(cbind(0:16, c("water","evergreen needleleaf forest","evergreen broadleaf forest","deciduous needleleaf forest","deciduous broadleaf forest","mixed forests","closed shrubland","open shrublands","woody savannas","savannas","grasslands","permanent wetlands","croplands","urban and built-up","cropland/natural vegetation mosaic","snow and ice","barren or sparsely vegetated")))
colnames(key) <- c('SYM','LC')

den <- as.data.frame(freq(raster(as.matrix(a))))
rownames(den) <- c("water","evergreen needleleaf forest","evergreen broadleaf forest","deciduous needleleaf forest","deciduous broadleaf forest","mixed forests","closed shrubland","open shrublands","woody savannas","savannas","grasslands","permanent wetlands","croplands","urban and built-up","cropland/natural vegetation mosaic","snow and ice","barren or sparsely vegetated")

key[key$SYM %in% den$value,3] <- den$count

dim(a)
b <- am6*a
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
#work on adding plus one column each time? may not be worth the coding hours




# elevation test
setwd("/Volumes/AOP-NEON1.4/VIP/")
setwd("~/Downloads/")

e <- read.csv("e10g")

install.packages("elevatr")
library(elevatr)

data(lake)
plot(get_elev_raster(lake, z = 2))
examp_df <- data.frame(x = runif(10, min = -73, max = -71), y = runif(10, min = 41,max = 45))
prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
plot(get_elev_raster(examp_df,prj = prj_dd,z =5))

install.packages("dismo")
library(dismo)

setwd("/Volumes/AOP-NEON1.4/VIP/JetVIP2/")

library(raster)
getData('ISO3')

#future problem of europe having so many countries for getting data:
# library(raster)
# misc = list()
# misc$countries = c("ZAF", "LSO", "SWZ", "ZWE", "MOZ", "NAM", "BWA")
# ctry_shps = do.call("bind", lapply(misc$countries, 
#                                    function(x) getData('GADM', country=x, level=0)))

library(dplyr)
library(purrr)


elevation <- getData("alt", country = "CAN")

#x <- terrain(elevation, opt = c("slope", "aspect"), unit = "degrees")
plot(elevation)

elevation1 <- getData("alt", country = "USA")
elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP2/USA2_msk_alt.grd`@extent
elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP2/USA3_msk_alt.grd`@extent
elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP2/USA4_msk_alt.grd`@extent

elevation_a <- elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP/USA1_msk_alt.grd`
elevation_b <- elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP/USA2_msk_alt.grd`
elevation_c <- elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP/USA3_msk_alt.grd`
elevation_d <- elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP/USA4_msk_alt.grd`
#plot(elevation_a, add = T) #but different color scheme, so need to change for plotting on same map
#plot(elevation_b, add = T)

elevation2 <- getData("alt", country = "MEX")

a <- merge(elevation,elevation_a,elevation_b,elevation_c,elevation_d,elevation2)

a@extent

am6 <- as.matrix(am6)
am6rotate <- raster(am6[nrow(am6):1,]) #need to flip rotate 3 on yaxis
extent(am6rotate) <- extent(-120,-94,20,70)

#crop extent
lonmin1 <- -120
lonmax1 <- -94
newext <- c(lonmin1, lonmax1, 20, 70)

am6_elevation <- crop(a,newext)

am6_elevation@ncols
am6_elevation@nrows

crs(am6_elevation)
plot(am6_elevation)

am6_elev_resampled <- resample(am6_elevation,am6rotate) #bilinear interpolation 
am6_elev_resampled@extent

b <- am6_elev_resampled * am6

b <- am6_elev_resampled * abs(am6)
summary(b)
as.matrix(b)

c <- as.numeric(unlist(b))
c <- c[complete.cases(c)]
boxplot(c, main = "N v S Significantly early/ late multiplied by elevation")

diff <- colMeans(m1)-colMeans(a2df)
d <- as.numeric(unlist(diff))
d[d == "NaN"] = "NA"  
d <- as.numeric(d)
cd <- data.frame(cbind(c,d))
cd_narm <- cd[complete.cases(cd),]
plot(cd_narm$d,cd_narm$c, xlab = "DOY",ylab = "Elevation")
legend("topright",
       c("colMeans(N)-colMeans(ALL)","colMeans(S)-colMeans(ALL)"),
       pch=c(1,1), # gives the legend appropriate symbols (lines)
       col=c("black","red"),
       cex = 0.75) # gives the legend lines the correct color and width

diff <- colMeans(m2)-colMeans(a2df)
d <- as.numeric(unlist(diff))
d[d == "NaN"] = "NA"  
d <- as.numeric(d)
cd <- data.frame(cbind(c,d))
cd_narm <- cd[complete.cases(cd),]
points(cd_narm$d,cd_narm$c, col = "red")

#length(seq(20,70,0.05)) #1001
#length(seq(-120,-94,0.05)) #521

spTransform(am6_elevation)
##############################################################################################################
## Koppen Geiger
install.packages("rgdal")
library(rgdal)
library(raster)
#t083<- raster('Beck_KG_V1_present_0p5.tif')

cc <- raster('KOPPEN_GEO_05.tif')
extent(cc)
e <- extent(-120,-94,20,70) 
cccropped <- crop(cc,e) #ncols 520 nrows 1000 
cccropped@data@values[cccropped@data@values==255] <- 0

freq(raster(as.matrix(cccropped)))

key <- as.data.frame(cbind(0:31, c("NA","Af","Am","Aw/As","BWh","BWk","BSh","BSk","Csa","Csb","Csc","Cwa","Cwb","Cwc","Cfa","Cfb","Cfc","Dsa","Dsb","Dsc","Dsd","Dwa","Dwb","Dwc","Dwd","Dfa","Dfb","Dfc","Dfd","ET","EF","?")))
colnames(key) <- c('SYM','KG')

den <- as.data.frame(freq(raster(as.matrix(cccropped))))

key[key$SYM %in% den$value,3] <- den$count

dim(cccropped)
b <- am6*cccropped
freq(raster(as.matrix(b)))

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
#work on adding plus one column each time? may not be worth the coding hours




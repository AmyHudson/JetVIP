# Land Cover, Koppen-Greiger climate classification, and Elevation test

# I've subsetted the MODIS LC into AM, JA, and JF regions in matlab, saved as .mat files, and now read them into R
setwd("/Volumes/AOP-NEON1.4/VIP/")

#Import the results of the significance file -1 and 1

am6 <- read.csv("am6_ns_sig.csv")
am6[,1] <- NULL
library(R.matlab)
a <- readMat("MODISLC_AM_6.mat")# proj=GCTP_GEO
a <- a$land2
den <- freq(raster(as.matrix(a)))

dim(a)
b <- am6*a
freq(raster(as.matrix(b)))

num_pos <- freq(raster(as.matrix(b)))[which(freq(raster(as.matrix(b)))[,1]>0),]
num_neg <- freq(raster(as.matrix(b)))[which(freq(raster(as.matrix(b)))[,1]<0),]
# match the value column of the numerator with the value of the denominator, divide count of that num by count of denom
lc <- matrix(NA, ncol = 2, nrow = 17)

for (i in 1:17){
  if (den[i,1] == num_pos[i,1]){
    lc[i,1] <- den[i,1]
    lc[i,2] <- num_pos[i,2]/den[i,2]
  } 
}

#how to fill in the data frame with missing sequential values? will that help me then match/merge them and divide by the count?
num_neg[1,2]/den[17,2]


# perhaps should write as a percentage of the total count in the region <- a loop for 


#merge these stats into one table
key <- cbind(0:16, c("water","evergreen needleleaf forest","evergreen broadleaf forest","deciduous needleleaf forest","deciduous broadleaf forest","mixed forests","closed shrubland","open shrublands","woody savannas","savannas","grasslands","permanent wetlands","croplands","urban and built-up","cropland/natural vegetation mosaic","snow and ice","barren or sparsely vegetated"))
as.numeric(key[,1])

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
elevation <- getData("alt", country = "CAN")
#x <- terrain(elevation, opt = c("slope", "aspect"), unit = "degrees")
plot(elevation)

elevation1 <- getData("alt", country = "USA")
elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP2/USA2_msk_alt.grd`@extent
elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP2/USA3_msk_alt.grd`@extent
elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP2/USA4_msk_alt.grd`@extent

elevation_a <- elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP2/USA1_msk_alt.grd`
elevation_b <- elevation1$`/Volumes/AOP-NEON1.4/VIP/JetVIP2/USA2_msk_alt.grd`
#plot(elevation_a, add = T) #but different color scheme, so need to change for plotting on same map
#plot(elevation_b, add = T)

elevation2 <- getData("alt", country = "MEX")

a <- merge(elevation,elevation_a,elevation_b,elevation2)
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

## Koppen Geiger

library(raster)
#t083<- raster('Beck_KG_V1_present_0p5.tif')

cc <- raster('KOPPEN_GEO_05.tif')
extent(cc)
e <- extent(-120,-94,20,70) 
cccropped <- crop(cc,e) #ncols 520 nrows 1000 

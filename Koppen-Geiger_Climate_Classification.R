#tiff classification 
setwd("~/Documents/JetVIP/Beck_KG_V1")
library(raster)
#t083<- raster('Beck_KG_V1_present_0p5.tif')
t083<- raster('KOPPEN_GEO_05.tif')

plot(t083)
#crop to 20 to 60N

extent(t083)
e <- extent(-120,-94,20,70)
t083c <- crop(t083,e)
plot(t083c)
lim_prec <- c(4,5,6,7,11,12,13,21,22,24)
lim_temp <- c(8,9,10,14,15,16,17,18,19,20,23,24,25,26,27,28,29,30)
lim_rad <- c(1,2,3)
t <- t083c
#FFFFFF #white
# base map for composites
t@legend@colortable[1] <- "#FFFFFF"
t@legend@colortable[2:31] <- "#D3D3D3"


#temp limited
t@legend@colortable[lim_prec] <- "#FFFFFF"
t@legend@colortable[lim_temp] <- "#000000"
t@legend@colortable[lim_rad] <- "#FFFFFF"
t@legend@colortable[31] <- "#FFFFFF"

#t[t == lim_prec] <- NA #@nodatavalue
#t[t == lim_rad] <- NA

plot(t, legend.only = T)

t <- t083c
t@legend@colortable[lim_prec] <- "#000000"
t@legend@colortable[lim_temp] <- "#FFFFFF"
t@legend@colortable[lim_rad] <- "#FFFFFF"
t@legend@colortable[31] <- "#FFFFFF"

plot(t, legend.only = T)

#an annual classification system
#subannual? What is limiting to cold

#group types together?

summary(t083c)

#################




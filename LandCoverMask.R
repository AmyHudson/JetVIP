#MCD12C1.A2001001.006.2018053185512.hdf

library(ncdf4)
library(rgdal)
library(gdalUtils)
gdalinfo('MCD12C1.A2001001.006.2018053185512.hdf')
b <- nc_open('MCD12C1.A2001001.006.2018053185512.hdf')
library(raster)
b <- raster('MCD12C1.A2001001.006.2018053185512.hdf')


library(maptools) 
getKMLcoordinates(textConnection(system("unzip -p User/amyhudson/Documents/JetVIP/Majority_Land_Cover_Type_1_in_MCD12C1.A2001001.006.20180531855.kml", intern = TRUE))) 
readOGR('Majority_Land_Cover_Type_1_in_MCD12C1.A2001001.006.20180531855.kml')
ogrListLayers('Majority_Land_Cover_Type_1_in_MCD12C1.A2001001.006.20180531855.kmz')


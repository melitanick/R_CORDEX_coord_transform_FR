###########################################################
# Task 2.3
# Read nc CORDEX files  
# Use file already reprojected from rotated projection
###########################################################

###########################################################
# DHMZ example
# M. Percec Tadic
# April 4, 2017.
###########################################################

library(ncdf4)
library(ncdf4.helpers)
library(raster)
library(maptools)
library(rgdal)

###########################################################
#
# Control grid - French, 9000 m,  ascii
# used for defining the extent of the area
#
###########################################################

control.dir <- "./controldir/" #dir
control.f <- "frgridascci1.asc" #file
control.in <- paste0(control.dir, control.f) #full path

control.r <- raster(control.in) #raster
res(control.r); proj4string(control.r)

###########################################################
# Set desired projection and domain extent for
# projected grids - pr based on control grids
###########################################################

pr.proj4 <- CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 
pr.ext <- extent(control.r)
proj4string(control.r) <- pr.proj4
spplot(control.r, scales=list(draw=T, y = list(rot = 90), cex=0.8))

###########################################################
#
# input nc CORDEX file
#
###########################################################

in.dir <- "./idir/"
in.f <- "tas_FR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1a_3hr_197001010000-197012312100_lonlat_timestep1.nc"
in.in <- paste0(in.dir, in.f)

nc.l <- nc_open(in.in)
nc.get.variable.list(nc.l)
nc.get.proj4.string(nc.l, "tas")
# so we need to know it from ...

# TODO here we could execute CDO from R
#
r1 <- brick(in.in, var='tas')[[1]] # brick is faster than stack, but only from single file

#For CORDEX grid the latlon on sphere is specified by the model settings. 
#But it seems there is an issue with reprojecting. This needs to be checked.
spf.proj4 <- CRS("+proj=longlat +a=6371229. +b=6371229. +towgs84=0,0,0 +no_defs") 
wgs.proj4 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4.l <- list(spf.proj4, wgs.proj4)

proj4string(r1) <- proj4.l[[1]] # chose sphere
r2 <- r1
proj4string(r2) <- proj4.l[[2]] # chose WGS

res(r1) #0.11°x0.11°
res(r2)

###########################################################
# Reproject and write ascii
###########################################################

out.dir <- "./odir/"
#dir.create(out.dir)
pr1 <- projectRaster(from=r1, to=control.r, method='ngb',
       overwrite=TRUE, format='ascii', filename= paste0(out.dir, "spf.", in.f, names(r1),".ascii"))

pr2 <- projectRaster(from=r2, to=control.r, method='ngb',
       overwrite=TRUE, format='ascii', filename= paste0(out.dir, "wgs.", in.f, names(r2), ".ascii"))

###########################################################
# Administrative borders
###########################################################

Fra <- getData('GADM', country='FRA', level=1)
class(Fra); proj4string(Fra)
Fra.pr <- spTransform(Fra, pr.proj4)
Fra.l <- as(Fra.pr, "SpatialLines") 

l1 <- list("sp.lines", Fra.l, col="black")
plot(Fra.l)

###########################################################
# Control points - French, 9000 m, shp
###########################################################

s1 <- readShapePoints(paste0(control.dir, "fr_natint1.shp"))
proj4string(s1) <- pr.proj4
l2 <- list("sp.points", s1, col="black", pch = 20, cex = 0.2)

###########################################################
# spplot
###########################################################

spplot(pr1, sp.layout=list(l1, l2), scales=list(draw=T, y = list(rot = 90), cex=0.8), main=names(r1))
spplot(pr2, sp.layout=list(l1, l2), scales=list(draw=T, y = list(rot = 90), cex=0.8), main=names(r2))

###########################################################
# Compare with files from IG
###########################################################

compare.dir <- "../coord_transform_IG/"
compare.f <-  "test2.tif"
compare.in <- paste0(compare.dir, compare.f)

compare.r <- brick(compare.in)

compare.b <- brick(pr1, pr2, compare.r)

spplot(compare.b[[c(1,3)]], sp.layout = list(l1, l2), scales=list(draw=T, y = list(rot = 90), cex=0.8), main="R  vs sh")

spplot(compare.b[[c(2,3)]], sp.layout = list(l1, l2), scales=list(draw=T, y = list(rot = 90), cex=0.8), main="R vs sh")

res(compare.b) #9000

save.image("06_crdtrans_CORDEX2EUCIRCLEnc.rda")

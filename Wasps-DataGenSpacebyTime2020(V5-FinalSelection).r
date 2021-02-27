### SCRIPT GENERATING DATA FOR THE OPTIMAL DESIGN ALGORITHM
##### Data: http://maps.elie.ucl.ac.be/CCI/viewer/index.php ### Water bodies 4.0

### Webpage https://mgimond.github.io/Spatial/point-pattern-analysis-in-r.html

##### Data: http://maps.elie.ucl.ac.be/CCI/viewer/index.php ### Water bodies 4.0

### Webpage https://mgimond.github.io/Spatial/point-pattern-analysis-in-r.html

##### Load the libraries
library(raster)
library(geosphere)
library(maptools)
library(MASS)
library(rgdal)
library(SDMTools)
library(spatstat)
library(sf)
library(beepr)
library(SDraw)
library(rgeos)


#### Import grid to use as a window
shape.grid<-readOGR("e:/CONTAIN/Wasps/Wasps2020/Grid_sel_final.shp")

shape.grid

plot(shape.grid)

### Centroids of each cell
cent.grid<-coordinates(shape.grid)

#### Location of CEHUm
x.cehum<-649573.53

y.cehum<-5593520.26

points(x.cehum, y.cehum, pch=19)

#### Euclidean distance between CEHUM and each cell
distance.tot.cehum<-round(apply(cent.grid, 1, function(x) sqrt((x[1]-x.cehum)^2+(x[2]-y.cehum)^2)), digits=0)

### Define the geographical window from the grid
win.grid<-as.owin(as.vector(extent(shape.grid)))

###### Import wasp data
shape.wasps<-readOGR("e:/CONTAIN/Wasps/Wasps2020/WaspColonies2020.shp")

shape.wasps

#### Altogether

plot(shape.grid)

points(shape.wasps, pch=19, col="yellow")

### Transform to point process

wasp.ppp<-ppp(coordinates(shape.wasps)[, 1], coordinates(shape.wasps)[, 2],  window=win.grid)

wasp.ppp

### 1st nearest neighbour

median(nndist(wasp.ppp, k=1))

### Create raster
res.raster<-1000          ### Resolution - 1000 metres in length

ext.ras<-extent(shape.grid)+1000        ### Plus a buffer of 1000-metres

### Number of columns of the raster (x)
ras.col<-round((ext.ras[2]-ext.ras[1])/res.raster, digits=0)

ras.col

### Number of rows of the raster (y)
ras.row<-round((ext.ras[4]-ext.ras[3])/res.raster, digits=0)

ras.row

### Create empty raster
ras.wasps<-raster(ext=ext.ras, nrow=ras.row, ncol=ras.col)

crs(ras.wasps)<-"+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

values(ras.wasps)<-1

### Check a plot

plot(ras.wasps)

points(shape.wasps, pch=19, col="black")

#### Rasterize
r.pres<-rasterize(shape.wasps, ras.wasps, field=rep(1, length(shape.wasps$ID_NEST)), fun='count', background=0)

values(r.pres)

summary(values(r.pres))

plot(r.pres)

plot(shape.wasps, add=TRUE)

########### TRACKS EFFORT
### Tracks
effort.wasps<-readOGR("e:/CONTAIN/Wasps/Wasps2020/Tracks.shp")

effort.wasps

plot(effort.wasps)

### Subset by survey type
random.search<-effort.wasps[effort.wasps$SEARCH=="Random", ]

semi.search<-effort.wasps[effort.wasps$SEARCH=="Semi-targeted", ]

target.search<-effort.wasps[effort.wasps$SEARCH=="Targeted", ]

### Rasters
rsp<-rasterToPolygons(ras.wasps)
crs(rsp)<-"+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

### Random search
rast.random<-intersect(random.search, rsp)

plot(rast.random)

rast.random$transect.length<-gLength(rast.random, byid=TRUE)

effort.random<-rasterize(rast.random, ras.wasps, field=rep(1, length(rast.random$transect.length)), fun='sum', background=0)

summary(effort.random)

plot(effort.random)

### Semi-targeted search
rast.semi<-intersect(semi.search, rsp)

plot(rast.semi)

rast.semi$transect.length<-gLength(rast.semi, byid=TRUE)

effort.semi<-rasterize(rast.semi, ras.wasps, field=rep(1, length(rast.semi$transect.length)), fun='sum', background=0)

summary(effort.semi)

plot(effort.semi)


### Targeted search
rast.target<-intersect(target.search, rsp)

plot(rast.target)

rast.target$transect.length<-gLength(rast.target, byid=TRUE)

effort.target<-rasterize(rast.target, ras.wasps, field=rep(1, length(rast.target$transect.length)), fun='sum', background=0)

summary(effort.target)

plot(effort.target)

####### Extract values to data frame
det.ab<-data.frame(n.nest=values(r.pres), random.eff=values(effort.random), semi.eff=values(effort.semi), target.eff=values(effort.target),
        tot.eff=values(effort.random)+values(effort.semi)+values(effort.target), dist.cehum=distance.tot.cehum)

######################################################## COVARIATES
#### Water stuff - length of water in each cell
water.shape<-readOGR("e:/CONTAIN/Wasps/Wasps2020/cursos agua region (hidroregion)_SAG.shp")

water.shape

plot(water.shape)

#### Rasterize
r.water<-rasterize(water.shape, ras.wasps, field=water.shape$LENGTH, fun='sum', background=0)

values(r.water)

summary(values(r.water))

plot(r.water)

plot(water.shape, add=TRUE)

#### Export raster
writeRaster(r.water, "e:/CONTAIN/Wasps/pop density/WaterBody-Length.grd", format="raster", overwrite=TRUE)

#### Add values to data frame
det.ab$lenght.water<-values(r.water)

########################### Population density 2020: from https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
pop.dens<-raster("e:/CONTAIN/Wasps/pop density/gpw_v4_population_density_rev11_2020_30_sec_6.asc")

pop.dens<-projectRaster(pop.dens, crs="+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" )

plot(pop.dens)

pop.dens

new.dens<-resample(pop.dens, ras.wasps, "bilinear")

dens2<-crop(pop.dens, extent(ras.wasps))

new.dens<-mask(new.dens, ras.wasps)

plot(new.dens)

writeRaster(new.dens, "e:/CONTAIN/Wasps/pop density/WaspsPopDens.grd", format="raster", overwrite=TRUE)

#### Add values to data frame
det.ab$pop.dens<-values(new.dens)

##############################################   LAND USES - FROM MAGDA
land.cov<-readOGR("e:/CONTAIN/Wasps/Wasps2020/usos.shp")

### 8 types of land use
use.class<-levels(land.cov$Uso_final)

use.class

#### Rasterize
proj.cov<-rasterize(land.cov, ras.wasps, field=land.cov$Uso_final, fun='max', background=NA)

values(proj.cov)

summary(values(proj.cov))

plot(proj.cov)

plot(proj.cov, add=TRUE)

#### Export raster
writeRaster(proj.cov, "e:/CONTAIN/Wasps/pop density/Land-Use.grd", format="raster", overwrite=TRUE)

#### Add values to data frame
det.ab$land.use<-use.class[values(proj.cov)]


################################################ NDVI https://land.copernicus.eu/global/products/ndvi
x.min<--74.111

x.max<--71.169

y.min<--41.758

y.max<--38.104

chile.ext<-extent(x.min, x.max, y.min, y.max)

#### SEPTEMBER
ndvi1<-raster("e:/CONTAIN/Wasps/pop density/NDVI-Sept2019.nc")

ndvi1<-crop(ndvi1, chile.ext)

ndvi1

plot(ndvi1)

proj.cover<-projectRaster(ndvi1, crs="+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" )

plot(proj.cover)

proj.cover

new.cover<-resample(proj.cover, ras.wasps, "bilinear")

cover2<-crop(proj.cover, extent(ras.wasps))

new.ndvi1<-mask(new.cover, ras.wasps)

plot(new.ndvi1)

writeRaster(new.ndvi1, "e:/CONTAIN/Wasps/pop density/NDVI-Sept19LosRios.grd", format="raster", overwrite=TRUE)

#### OCTOBER
ndvi2<-raster("e:/CONTAIN/Wasps/pop density/NDVI-Oct2019.nc")

ndvi2<-crop(ndvi2, chile.ext)

ndvi2

plot(ndvi2)

proj.cover<-projectRaster(ndvi2, crs="+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" )

plot(proj.cover)

proj.cover

new.ndvi2<-resample(proj.cover, ras.wasps, "bilinear")

cover2<-crop(proj.cover, extent(ras.wasps))

new.ndvi2<-mask(new.ndvi2, ras.wasps)

plot(new.ndvi2)

writeRaster(new.ndvi2, "e:/CONTAIN/Wasps/pop density/NDVI-Oct19LosRios.grd", format="raster", overwrite=TRUE)


#### NOVEMBER
ndvi3<-raster("e:/CONTAIN/Wasps/pop density/NDVI-Nov2019.nc")

ndvi3<-crop(ndvi3, chile.ext)

ndvi3

plot(ndvi3)

proj.cover<-projectRaster(ndvi3, crs="+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" )

plot(proj.cover)

proj.cover

new.ndvi3<-resample(proj.cover, ras.wasps, "bilinear")

cover3<-crop(proj.cover, extent(ras.wasps))

new.ndvi3<-mask(new.ndvi3, ras.wasps)

plot(new.ndvi3)

writeRaster(new.ndvi3, "e:/CONTAIN/Wasps/pop density/NDVI-Nov19LosRios.grd", format="raster", overwrite=TRUE)

#### DCEMBER
ndvi4<-raster("e:/CONTAIN/Wasps/pop density/NDVI-Dec2019.nc")

ndvi4<-crop(ndvi4, chile.ext)

ndvi4

plot(ndvi4)

proj.cover<-projectRaster(ndvi4, crs="+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" )

plot(proj.cover)

proj.cover

new.ndvi4<-resample(proj.cover, ras.wasps, "bilinear")

cover4<-crop(proj.cover, extent(ras.wasps))

new.ndvi4<-mask(new.ndvi4, ras.wasps)

plot(new.ndvi4)

writeRaster(new.ndvi4, "e:/CONTAIN/Wasps/pop density/NDVI-Dec19LosRios.grd", format="raster", overwrite=TRUE)

#### JANUARY 2020
ndvi5<-raster("e:/CONTAIN/Wasps/pop density/NDVI-Jan20.nc")

ndvi5<-crop(ndvi5, chile.ext)

ndvi5

plot(ndvi5)

proj.cover<-projectRaster(ndvi5, crs="+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" )

plot(proj.cover)

proj.cover

new.ndvi5<-resample(proj.cover, ras.wasps, "bilinear")

cover5<-crop(proj.cover, extent(ras.wasps))

new.ndvi5<-mask(new.ndvi5, ras.wasps)

plot(new.ndvi5)

writeRaster(new.ndvi5, "e:/CONTAIN/Wasps/pop density/NDVI-Jan20LosRios.grd", format="raster", overwrite=TRUE)

#### Add values to data frame
ndvi.tot<-data.frame(values(new.ndvi1), values(new.ndvi2), values(new.ndvi3), values(new.ndvi4), values(new.ndvi5))

det.ab$mean.ndvi<-rowMeans(ndvi.tot)

det.ab$median.ndvi<-apply(ndvi.tot, 1, median)

det.ab$var.ndvi<-apply(ndvi.tot, 1, var)

det.ab

summary(det.ab)

### Export

write.table(det.ab, "e:/CONTAIN/Wasps/Wasps2020/Nest-Count-OptimalExperimental.csv", sep=",")

write.table(det.ab, "e:/CONTAIN/Experimental Design/Wasps/Nest-Count-OptimalExperimental.csv", sep=",")

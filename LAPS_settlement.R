# Load libraries
library(sf)
library(sp)
library(RANN)
library(dplyr)
library(raster)
library(ncdf4)

#-------------------------------------
# Make Settlement Files
#-------------------------------------
# Read in the habitat file
landxy <- read.csv('HIreefsNew.csv', header=F)
names(landxy) <- c('lat', 'lon', 'cover', 'island')
# Remove points that have less than 1% coral cover
landxy <- filter(landxy, cover >= 0.01, lat > 18)

# Ocean points
# Make an evenly spaced grid
lon <- seq(120,240,by=0.01)
lat <- seq(0,50,by=0.01)
lonlat <- expand.grid(lon=lon, lat=lat)

# find nearest habitat
# fast nearest neighbour search for any habitat pixel within 0.05 degrees of each grid point
# Points outside the search radius will be labeled with a 0 and a length of 1.340781e+154 degrees
closest <- nn2(landxy[,2:1], lonlat, k=1, searchtype = "radius", radius=0.05)
# If you want all distances to nearest settlement habitat, just leave the searchtype and radius out like below
# closest <- nn2(landxy[,2:1], lonlat, k=1)
closest <- sapply(closest, cbind) %>% as_tibble

# add nearest neighbor info to ocean data
test <- cbind(lonlat, closest)
names(test) <- c('Longitude', 'Latitude', 'IslandIndex', 'Distance')
# Presence/absence of settlement habitat
test$Settle <- ifelse(test$IslandIndex == 0, 0, 1)
# Remove the arbitrary distance nn2 adds to points outside the search radius (only relevant if you did the radius search on 'closest')
test$Distance <- ifelse(test$Distance >= 1.3407e154, NA, test$Distance)

# make into raster for easier saving
# I neded up making individual rastes and saving as indivudual ncdf files since the layers didn't do well during saving and importing into Python
testRas <- rasterFromXYZ(test[,c(1,2,5)], crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
testRasIsl <- rasterFromXYZ(test[,c(1,2,3)], crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
testRasDist <- rasterFromXYZ(test[,c(1,2,4)], crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# Write files
writeRaster(testRas, 'LAPS_settlement_habitat_Settle.nc', overwrite=T)
writeRaster(testRasIsl, 'LAPS_settlement_habitat_IslIndex.nc', overwrite=T)
writeRaster(testRasDist, 'LAPS_settlement_habitat_Dist.nc', overwrite=T)


#-------------------------------------
# Check Parcels output
#-------------------------------------
# Read trajectory file
dat <- nc_open('LAPS_sensitivity_n1_pld180_20m_Kh0_nday545_daily_07012012.nc')
lo <- ncvar_get(dat, 'lon')
la <- ncvar_get(dat, 'lat')
a <- ncvar_get(dat, 'age')
i <- ncvar_get(dat, 'trajectory')
site <- ncvar_get(dat, 'releaseSite')
t  <- ncvar_get(dat, 'time')
tunit <- dat$var$time$units
d <- ncvar_get(dat, 'distance')
# If you have the settle variable read it in here, otherwise skip it (obviously)!
#s <- ncvar_get(dat, 'settle')
nc_close(dat)

# Convert into vectors and merge data
lon <- as.vector(lo)
lat <- as.vector(la)
tvec <- as.vector(t)
age <- as.vector(a)
id <- as.vector(i)
site <- as.vector(site)
distance <- as.vector(d)
#settle <- as.vector(s)

# Fix dates
tustr <- strsplit(tunit, " ")
stime <- as.POSIXct(unlist(tustr)[3])
time <- as.POSIXct(tvec, origin=stime)

# Make trajectory dataframe
traj <- data.frame(id, lon, lat, age, time, site, distance)#, settle)

#-------------------------------------
# What particles settled and where?
# Extract the value in the settlement raster we made above underneath each saved trajectory point
# You can do this with distance to settlement habitat raster (as shown below) or if you are only interested in habitat pixel use the island index raster
traj$settleDistance <- extract(testRasDist, traj[,2:3])
# Save only settled particles that are older than the precompetency period (in this example that's 14 days) and are within a certain distance from the habitat (in degrees)
trajSettle <- traj %>%
  filter(settleDistance >= 0.05, age >= (86400*14))
# Then you can plot the trajectories from the settled particles, calculate mean distance traveled for settled particles, see what sites contributed to settlement

#-------------------------------------
# If you have a settle column in your Parcels output, to find the site that belongs to it do this:
trajSettle <- traj %>% 
  filter(settle !=0, age >= (86400*14))
traj$settleHabPixel <- extract(testRasIsl, trajSettle[,1:2])

#-------------------------------------
# Or, you can just directly calculate the distance between the trajectory points and the nearest habitat skip the rasters and do this instead
traj <- na.omit(traj)
closest2 <- nn2(landxy[,2:1], traj[,2:3], k=1, searchtype = "radius", radius=0.05)
closest2 <- sapply(closest2, cbind) %>% 
  as_tibble %>% 
  rename(IslandIndex=nn.idx, habDistance=nn.dists) %>% 
  bind_cols(traj) %>% 
  dplyr::select(id:distance, IslandIndex, habDistance)



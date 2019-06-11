# Load and save raw data internally to package simfish
# Do not run this as a script


#install.packages(c("rgeos", "maptools", "foreign")) # dependencies
library(tidyverse)
library(sp)
library(PBSmapping)
library(rgdal)
library(raster)
library(automap)
library(gstat)
library(marmap)
library(maptools)



# Do once for the package
# devtools::use_data_raw()

# Load csv with raw data for the simulation
All_data <- read.csv(paste0(getwd(), "/data-raw/catch0.csv"))

# configure according to Kotaro's original script
All_data$CPUE <- with(All_data, WEIGHT/(DISTANCE_FISHED * NET_WIDTH * .001)) # CPUE in kg/km^2
All_data$LONGITUDE = All_data$START_LONGITUDE
All_data$LATITUDE = All_data$START_LATITUDE

#table(All_data$SPECIES_CODE)

# Change projection to UTM coordinates to make sure that distance are preserved betwen locations approximately
Data_LL <- All_data[,c('LONGITUDE', 'LATITUDE')]
coordinates(Data_LL) <- c("LONGITUDE", "LATITUDE")
proj4string(Data_LL) <- CRS("+proj=longlat +datum=WGS84")
Data_UTM <- spTransform(Data_LL, CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

All_data <- data.frame(All_data[,c('YEAR', 'SURVEY_NAME', 'SPECIES_CODE', 'BOTTOM_DEPTH', 'SURFACE_TEMPERATURE', 'GEAR_TEMPERATURE', 'CPUE')], Data_UTM@coords)
All_data$LONGITUDE <- All_data$LONGITUDE
All_data$LATITUDE <- All_data$LATITUDE


### Prepare the prediction grid
### WARNING: Prediction/extrapolation of environmental covariate (i.e. temperature) is based on the survey data.
### The approach is based on fitting the best variogram model for each covariate and year, then extrapolate to other areas
### This means that for years without survey, these predictions could be pretty vague.

# Narrow prediction extent to nnmfs area polygon shapefile
simplenmfs <- readOGR(dsn="data/AKmaps_bathy",layer="simplenmfs") 

# get bathymetry for prediction extent, working around antimeridian
bathy1 <- getNOAA.bathy(lon1 = extent(simplenmfs)@xmin, lon2 = 180, lat1 = extent(simplenmfs)@ymin, lat2 = extent(simplenmfs)@ymax, resolution = 2)
bathy2 <- getNOAA.bathy(lon1 = -180, lon2 = -130, lat1 = extent(simplenmfs)@ymin, lat2 = extent(simplenmfs)@ymax, resolution = 2)
temp <- fortify.bathy(bathy2)
temp$x[temp$x<0] <- temp$x[temp$x<0] + 360
Predict_data <- bind_rows(fortify.bathy(bathy1),temp)
Predict_data$x <- Predict_data$x-360
Predict_data = as.raster(as.bathy(Predict_data))
Predict_data = projectRaster(Predict_data, crs = "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
#> res(Predict_data)
#[1] 1940 3690 # 2 arc minute
#[1] 1460 2770 # 1.5 arc minute
Predict_data = as.xyz(as.bathy(Predict_data))
colnames(Predict_data) <- c("X", "Y", "BOTTOM_DEPTH")
Predict_data$BOTTOM_DEPTH <- -Predict_data$BOTTOM_DEPTH


    ### Now import shapefiles from different AK region, reproject to same projection method, then use it to define prediction for each specific regions.
    # EBS shelf
    EBS_shelf <- readOGR(dsn=paste0(getwd(), "/data/shapefiles/EBS_shelfStrataNoland.shp"))
    EBS_shelf <- spTransform(EBS_shelf, CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
    lps <- getSpPPolygonsLabptSlots(EBS_shelf)
    IDOneBin <- cut(lps[,1], range(lps[,1]), include.lowest=TRUE)
    NcDissolve   <- unionSpatialPolygons(EBS_shelf ,IDOneBin)
    EBS_shelfOne <- SpatialPolygons2PolySet(NcDissolve)

    # EBS slope
    EBS_slope <- readOGR(dsn=paste0(getwd(), "/data/shapefiles/EBS_slopeStrata.shp"))
    EBS_slope <- spTransform(EBS_slope, CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
    lps <- getSpPPolygonsLabptSlots(EBS_slope)
    IDOneBin <- cut(lps[,1], range(lps[,1]), include.lowest=TRUE)
    NcDissolve   <- unionSpatialPolygons(EBS_slope ,IDOneBin)
    EBS_slopeOne <- SpatialPolygons2PolySet(NcDissolve)

    # BS north
    BS_north <- readOGR(dsn=paste0(getwd(), "/data/shapefiles/NBS_shelfStrataNoland.shp"))
    BS_north <- spTransform(BS_north, CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
    lps <- getSpPPolygonsLabptSlots(BS_north)
    IDOneBin <- cut(lps[,1], range(lps[,1]), include.lowest=TRUE)
    NcDissolve   <- unionSpatialPolygons(BS_north ,IDOneBin)
    BS_northOne <- SpatialPolygons2PolySet(NcDissolve)

    # GOA
		GOA <- readOGR(dsn=paste0(getwd(), "/data/shapefiles/GOA_erase.shp"))
    GOA <- spTransform(GOA, CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
    GOAOne <- SpatialPolygons2PolySet(GOA)

		# AI
		AI <- readOGR(dsn=paste0(getwd(), "/data/shapefiles/AI_dissolved_noland.shp"))
    AI <- spTransform(AI, CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
    AIOne <- SpatialPolygons2PolySet(AI)

    # Merging EBS shelf and slope
    EBS <- union(EBS_shelf, EBS_slope)
    lps <- getSpPPolygonsLabptSlots(EBS)
    IDOneBin <- cut(lps[,1], range(lps[,1]), include.lowest=TRUE)
    NcDissolve   <- unionSpatialPolygons(EBS ,IDOneBin)
    EBSOne <- SpatialPolygons2PolySet(NcDissolve)

    # Merging EBS shelf and north
    EBS_shelf_north <- union(EBS_shelf, BS_north)
    lps <- getSpPPolygonsLabptSlots(EBS_shelf_north)
    IDOneBin <- cut(lps[,1], range(lps[,1]), include.lowest=TRUE)
    NcDissolve   <- unionSpatialPolygons(EBS_shelf_north ,IDOneBin)
    EBS_shelf_northOne <- SpatialPolygons2PolySet(NcDissolve)


    ### Now include the SST values based on surface kriging of available temperature data (BS+GOA)
    dat_new <- subset(All_data, select=c("YEAR", "SURFACE_TEMPERATURE", "GEAR_TEMPERATURE", "LATITUDE", "LONGITUDE"))
    dat_new <- dat_new[!duplicated(dat_new),]
    coordinates(dat_new) <- ~ LONGITUDE + LATITUDE

    YEARS = 1982:2016
    raster_temp_surf <- list()
    raster_temp_bottom <- list()
    for (yr in seq_along(YEARS))
    {
      # For surface temperature
      dat1 <- subset(dat_new, subset=c(YEAR == YEARS[yr]))
      dat1 <- dat1[!is.na(dat1$SURFACE_TEMPERATURE),]
      m <- autofitVariogram(SURFACE_TEMPERATURE~1, dat1)
      plot(m)
      v <- variogram(SURFACE_TEMPERATURE~1, dat1)	#create a variogram of the sorting data
      m <- fit.variogram(v, vgm(psill=m$var_model[2,2], model=as.character(m$var_model[2,1]), range=m$var_model[2,3], nugget =m$var_model[1,2], kappa=m$var_model[2,4])) #fit a model to the variogram
      # GETTING SOME CONVERGENCE WARNINGS FROM FIT.VARIOGRAM ABOVE, BUT PLOTTED VARIOGRAMS LOOK DECENT
      plot(v, model= m)

      bathy <- rasterFromXYZ(Predict_data[,c(1,2)])
      g <- gstat(id = "SURFACE_TEMPERATURE", formula = SURFACE_TEMPERATURE~1, data=dat1, model = m, nmax=5)
      raster_temp_surf[[yr]]  <- raster::interpolate(bathy, g, xyOnly=TRUE, progress="text",overwrite=TRUE ) #Interpolate the object to a raster

      # For gear temperature
      dat1 <- subset(dat_new, subset=c(YEAR == YEARS[yr]))
      dat1 <- dat1[!is.na(dat1$GEAR_TEMPERATURE),]
      m <- autofitVariogram(GEAR_TEMPERATURE~1, dat1)
      plot(m)
      v <- variogram(GEAR_TEMPERATURE~1, dat1)	#create a variogram of the sorting data
      m <- fit.variogram(v, vgm(psill=m$var_model[2,2], model=as.character(m$var_model[2,1]), range=m$var_model[2,3], nugget=m$var_model[1,2], kappa=m$var_model[2,4]))    #fit a model to the variogram
      plot(v, model= m)

      bathy <- rasterFromXYZ(Predict_data[,c(1,2)])
      g <- gstat(id = "GEAR_TEMPERATURE", formula = GEAR_TEMPERATURE~1, data=dat1, model = m, nmax=5)
      raster_temp_bottom[[yr]] <- raster::interpolate(bathy, g, xyOnly=TRUE, progress="text", overwrite=TRUE) #Interpolate the object to a raster
    }

    bathy <- rasterFromXYZ(Predict_data)

    # save the data for the whole area
    Grids <- expand.grid(LONG=sort(unique(Predict_data$X)), LAT=sort(unique(Predict_data$Y)))
    Predict_data <- data.frame(Grids, BOTTOM_DEPTH = extract(bathy, Grids), SURFACE_TEMPERATURE1982=extract(raster_temp_surf[[1]], Grids), SURFACE_TEMPERATURE1983= extract(raster_temp_surf[[2]], Grids), SURFACE_TEMPERATURE1984= extract(raster_temp_surf[[3]], Grids), SURFACE_TEMPERATURE1985= extract(raster_temp_surf[[4]], Grids), SURFACE_TEMPERATURE1986= extract(raster_temp_surf[[5]], Grids), SURFACE_TEMPERATURE1987= extract(raster_temp_surf[[6]], Grids), SURFACE_TEMPERATURE1988= extract(raster_temp_surf[[7]], Grids), SURFACE_TEMPERATURE1989= extract(raster_temp_surf[[8]], Grids), SURFACE_TEMPERATURE1990= extract(raster_temp_surf[[9]], Grids), SURFACE_TEMPERATURE1991= extract(raster_temp_surf[[10]], Grids), SURFACE_TEMPERATURE1992= extract(raster_temp_surf[[11]], Grids), SURFACE_TEMPERATURE1993= extract(raster_temp_surf[[12]], Grids), SURFACE_TEMPERATURE1994= extract(raster_temp_surf[[13]], Grids), SURFACE_TEMPERATURE1995= extract(raster_temp_surf[[14]], Grids), SURFACE_TEMPERATURE1996= extract(raster_temp_surf[[15]], Grids), SURFACE_TEMPERATURE1997= extract(raster_temp_surf[[16]], Grids), SURFACE_TEMPERATURE1998= extract(raster_temp_surf[[17]], Grids), SURFACE_TEMPERATURE1999= extract(raster_temp_surf[[18]], Grids), SURFACE_TEMPERATURE2000= extract(raster_temp_surf[[19]], Grids), SURFACE_TEMPERATURE2001= extract(raster_temp_surf[[20]], Grids), SURFACE_TEMPERATURE2002= extract(raster_temp_surf[[21]], Grids), SURFACE_TEMPERATURE2003= extract(raster_temp_surf[[22]], Grids), SURFACE_TEMPERATURE2004= extract(raster_temp_surf[[23]], Grids), SURFACE_TEMPERATURE2005= extract(raster_temp_surf[[24]], Grids), SURFACE_TEMPERATURE2006= extract(raster_temp_surf[[25]], Grids), SURFACE_TEMPERATURE2007= extract(raster_temp_surf[[26]], Grids), SURFACE_TEMPERATURE2008= extract(raster_temp_surf[[27]], Grids), SURFACE_TEMPERATURE2009= extract(raster_temp_surf[[28]], Grids), SURFACE_TEMPERATURE2010= extract(raster_temp_surf[[29]], Grids), SURFACE_TEMPERATURE2011= extract(raster_temp_surf[[30]], Grids), SURFACE_TEMPERATURE2012= extract(raster_temp_surf[[31]], Grids), SURFACE_TEMPERATURE2013= extract(raster_temp_surf[[32]], Grids), SURFACE_TEMPERATURE2014= extract(raster_temp_surf[[33]], Grids), SURFACE_TEMPERATURE2015= extract(raster_temp_surf[[34]], Grids), SURFACE_TEMPERATURE2016= extract(raster_temp_surf[[35]], Grids),
                               GEAR_TEMPERATURE1982= extract(raster_temp_bottom[[1]], Grids), GEAR_TEMPERATURE1983= extract(raster_temp_bottom[[2]], Grids), GEAR_TEMPERATURE1984= extract(raster_temp_bottom[[3]], Grids), GEAR_TEMPERATURE1985= extract(raster_temp_bottom[[4]], Grids), GEAR_TEMPERATURE1986= extract(raster_temp_bottom[[5]], Grids), GEAR_TEMPERATURE1987= extract(raster_temp_bottom[[6]], Grids), GEAR_TEMPERATURE1988= extract(raster_temp_bottom[[7]], Grids), GEAR_TEMPERATURE1989= extract(raster_temp_bottom[[8]], Grids), GEAR_TEMPERATURE1990= extract(raster_temp_bottom[[9]], Grids), GEAR_TEMPERATURE1991= extract(raster_temp_bottom[[10]], Grids), GEAR_TEMPERATURE1992= extract(raster_temp_bottom[[11]], Grids), GEAR_TEMPERATURE1993= extract(raster_temp_bottom[[12]], Grids), GEAR_TEMPERATURE1994= extract(raster_temp_bottom[[13]], Grids), GEAR_TEMPERATURE1995= extract(raster_temp_bottom[[14]], Grids), GEAR_TEMPERATURE1996= extract(raster_temp_bottom[[15]], Grids), GEAR_TEMPERATURE1997= extract(raster_temp_bottom[[16]], Grids), GEAR_TEMPERATURE1998= extract(raster_temp_bottom[[17]], Grids), GEAR_TEMPERATURE1999= extract(raster_temp_bottom[[18]], Grids), GEAR_TEMPERATURE2000= extract(raster_temp_bottom[[19]], Grids), GEAR_TEMPERATURE2001= extract(raster_temp_bottom[[20]], Grids), GEAR_TEMPERATURE2002= extract(raster_temp_bottom[[21]], Grids), GEAR_TEMPERATURE2003= extract(raster_temp_bottom[[22]], Grids), GEAR_TEMPERATURE2004= extract(raster_temp_bottom[[23]], Grids), GEAR_TEMPERATURE2005= extract(raster_temp_bottom[[24]], Grids), GEAR_TEMPERATURE2006= extract(raster_temp_bottom[[25]], Grids), GEAR_TEMPERATURE2007= extract(raster_temp_bottom[[26]], Grids), GEAR_TEMPERATURE2008= extract(raster_temp_bottom[[27]], Grids), GEAR_TEMPERATURE2009= extract(raster_temp_bottom[[28]], Grids), GEAR_TEMPERATURE2010= extract(raster_temp_bottom[[29]], Grids), GEAR_TEMPERATURE2011= extract(raster_temp_bottom[[30]], Grids), GEAR_TEMPERATURE2012= extract(raster_temp_bottom[[31]], Grids), GEAR_TEMPERATURE2013= extract(raster_temp_bottom[[32]], Grids), GEAR_TEMPERATURE2014= extract(raster_temp_bottom[[33]], Grids), GEAR_TEMPERATURE2015= extract(raster_temp_bottom[[34]], Grids), GEAR_TEMPERATURE2016= extract(raster_temp_bottom[[35]], Grids))

    # Now include the above regions definition into the prediction data
    ## EBS shelf
			if(length(unique(EBS_shelfOne$SID))>1)
			{
				EBS_shelfRegion <- rep(0,nrow(Predict_data))
				for (i in unique(EBS_shelfOne$SID))
				{
					if (all(diff(subset(EBS_shelfOne, SID==i)$POS)>=0)) 
					{
						temp <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(EBS_shelfOne, SID==i)$X, subset(EBS_shelfOne, SID==i)$Y, mode.checked=FALSE)
						EBS_shelfRegion[which(temp==1)] <- 1
					}
					if (all(diff(subset(EBS_shelfOne, SID==i)$POS)<=0))
					{
						temp <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(EBS_shelfOne, SID==i)$X, subset(EBS_shelfOne, SID==i)$Y, mode.checked=FALSE)
						EBS_shelfRegion[which(temp==1)] <- 0
					}	
				}
			} else {
				EBS_shelfRegion <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(EBS_shelfOne, SID==1)$X, subset(EBS_shelfOne, SID==1)$Y, mode.checked=FALSE)
			}
    ## EBS slope
			if(length(unique(EBS_slopeOne$SID))>1)
			{
				EBS_slopeRegion <- rep(0,nrow(Predict_data))
				for (i in unique(EBS_slopeOne$SID))
				{
					if (all(diff(subset(EBS_slopeOne, SID==i)$POS)>=0)) 
					{
						temp <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(EBS_slopeOne, SID==i)$X, subset(EBS_slopeOne, SID==i)$Y, mode.checked=FALSE)
						EBS_slopeRegion[which(temp==1)] <- 1
					}
					if (all(diff(subset(EBS_slopeOne, SID==i)$POS)<=0))
					{
						temp <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(EBS_slopeOne, SID==i)$X, subset(EBS_slopeOne, SID==i)$Y, mode.checked=FALSE)
						EBS_slopeRegion[which(temp==1)] <- 0
					}	
			  }
			} else {
				EBS_slopeRegion <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(EBS_slopeOne, SID==1)$X, subset(EBS_slopeOne, SID==1)$Y, mode.checked=FALSE)
			}
    ## EBS north
			if(length(unique(BS_northOne$SID))>1)
			{
				BS_northRegion <- rep(0,nrow(Predict_data))
				for (i in unique(BS_northOne$SID))
				{
					if (all(diff(subset(BS_northOne, SID==i)$POS)>=0)) 
					{
						temp <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(BS_northOne, SID==i)$X, subset(BS_northOne, SID==i)$Y, mode.checked=FALSE)
						BS_northRegion[which(temp==1)] <- 1
					}
					if (all(diff(subset(BS_northOne, SID==i)$POS)<=0)) 
					{
						temp <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(BS_northOne, SID==i)$X, subset(BS_northOne, SID==i)$Y, mode.checked=FALSE)
						BS_northRegion[which(temp==1)] <- 0
					}		
				}
			} else {
				BS_northRegion <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(BS_northOne, SID==1)$X, subset(BS_northOne, SID==1)$Y, mode.checked=FALSE)
			}
    ## AI
    # AIRegion_new <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(AInewOne, SID==1)$X, subset(AInewOne, SID==1)$Y, mode.checked=FALSE)
			if(length(unique(AIOne$SID))>1)
			{
				AIRegion <- rep(0,nrow(Predict_data))
				for (i in unique(AIOne$SID))
				{
					if (all(diff(subset(AIOne, SID==i)$POS)>=0)) 
					{
						temp <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(AIOne, SID==i)$X, subset(AIOne, SID==i)$Y, mode.checked=FALSE)
						AIRegion[which(temp==1)] <- 1
					}
					if (all(diff(subset(AIOne, SID==i)$POS)<=0)) 
					{
						temp <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(AIOne, SID==i)$X, subset(AIOne, SID==i)$Y, mode.checked=FALSE)
						AIRegion[which(temp==1)] <- 0
					}		
				}
			} else {
				AIRegion <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(AIOne, SID==1)$X, subset(AIOne, SID==1)$Y, mode.checked=FALSE)
			}
    ### GOA
			if(length(unique(GOAOne$SID))>1)
			{
				GOARegion <- rep(0,nrow(Predict_data))
				for (i in unique(GOAOne$SID))
				{
					if (i==94) 	# this is the GOA survey area map without the land masses 
					{
					temp <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(GOAOne, SID==i)$X, subset(GOAOne, SID==i)$Y, mode.checked=FALSE)
					GOARegion[which(temp==1)] <- 1
					}
					if (i!=94) # this is all the land masses
					{
					temp <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(GOAOne, SID==i)$X, subset(GOAOne, SID==i)$Y, mode.checked=FALSE)
					GOARegion[which(temp==1)] <- 0
					}		
				}
			} else {
				GOARegion <- point.in.polygon(Predict_data$LONG, Predict_data$LAT, subset(GOAOne, SID==1)$X, subset(GOAOne, SID==1)$Y, mode.checked=FALSE)
			}

 #   Predict_data <- data.frame(Predict_data, EBS=EBS_shelfRegion, SLP=EBS_slopeRegion, NBS=BS_northRegion, GOA=GOARegion, AI=AIRegion, AI_old=AIRegion_new)
    Predict_data <- data.frame(Predict_data, EBS=EBS_shelfRegion, SLP=EBS_slopeRegion, NBS=BS_northRegion, GOA=GOARegion, AI=AIRegion)

    save(Predict_data, file=paste0(getwd(), "/data-raw/Predict_data.Rdata"))


devtools::use_data(All_data, Predict_data, internal = FALSE, overwrite = TRUE)

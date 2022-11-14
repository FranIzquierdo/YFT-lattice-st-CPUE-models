#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lattice data CPUE standardization NOAA YFT  #
# Prepare original dataset for INLA model     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Francisco Izquierdo #
#~~~~~~~~~~~~~~~~~~~~~~
# 30/10/2022 #        
#~~~~~~~~~~~~~

## press Ctrl + shift + O to see the script outline

# start here -------------------------------------------------------------------

rm(list=ls()) 

## read original dataset
df<-read.csv("./data/original_CPUE_grid_data.csv") 

## create outdir
out<-paste(getwd(),"/output", sep="")
dir.create(out)

## create dir
out_00<-paste(getwd(),"/output/00 prepare data", sep="")
dir.create(out_00)

# redefine vars ----------------------------------------------------------------

## define variable year and season from pseudoyear (81:256) column
## one year = 4 pseudo-years = 4 seasons
(nyears<-length(unique(df$year))/4) # 176 pseudoyear/4 seas = 44 years(1972:2015)

## create levels for year, season and combined yearseas
lyear<-rep(seq(from=1972, to=2015), each=4); length(lyear) # 176
lseas<-rep(1:4, times=nyears); length(lseas) # 176
lyearseas<-paste(lyear, lseas, sep=".")

## replicate factor columns in data frame
df$pseudoyear<-as.factor(df$year); levels(df$pseudoyear) 
df$year<-as.factor(df$year); levels(df$year)
df$seas<-as.factor(df$year); levels(df$seas)
df$yearseas<-as.factor(df$year); levels(df$yearseas)

## rename factor levels
levels(df$year) <- lyear
levels(df$seas)<-lseas
levels(df$yearseas)<-lyearseas

## prepare ID variables for inla (from 1:n)
df$regionID<-as.factor(df$region)
levels(df$regionID)<-c(1,1,2,3,4) ## region ID
df$yearID<-as.numeric(as.integer(df$year)) ## year ID
df$psyearID<-as.numeric(as.integer(df$pseudoyear))

# polygrids --------------------------------------------------------------------

## we have the coordinate centroids for each spatial cell polygon
## we need to construct cell polygons in order to fit a lattice spatial model

xy<-df[,7:6] # lat long
plot(xy)

## with this function we can create polygon grids from points
## note that dlon and dlat are degrees from the centroid of each point (5x5)

library(CCAMLRGIS)
MyGrid=create_PolyGrids(Input=df,dlon=5,dlat=5 , NamesIn=c('lat2','lon2'))
plot(st_geometry(MyGrid),main="WGS 84 / NSIDC EASE-Grid 2.0 South", axes=TRUE)

## set CRS = 4325 (WGS84)
MyGrid<-st_transform(MyGrid, 4326) # WGS84
png(paste0(out_00,"/00 poly gridsand centroid pts WGS84.png"), width = 1000, height = 800, res = 150)
plot(st_geometry(MyGrid),main="pts and poly grids (5x5 degrees) WGS84", axes=TRUE)
points(xy$lon2, xy$lat2,  col="blue")
dev.off()

# join  ------------------------------------------------------------------------

## we need to join both datasets, the spatial poly grid and the observed data
## the order of the grid cells must be the same between both datasets

## note that the function polygrids modifies the decimal values of the coordinates
## in order to do a join we need to have the same decimals, so we correct this:
MyGrid$Centrelon_r<-round(MyGrid$Centrelon, 1)
MyGrid$Centrelat_r<-ceiling(MyGrid$Centrelat)-0.5

## set same names of lat long variables to do the join
colnames(df)<-c("year","lat","lon","cpue","region","Centrelat_r","Centrelon_r",
                "sim","pseudoyear","seas", "yearseas", "regionID", "yearID","psyearID")
## join
library(dplyr)
dft<-left_join(MyGrid, df, by=c("Centrelat_r","Centrelon_r"))

## rename ID as polyID
colnames(dft)[colnames(dft) == 'ID'] <- 'polyID'
dft$polyID<-as.numeric(dft$polyID)

# poly grid --------------------------------------------------------------------

## save poly grid dataset as unique polyID and geometry 
poly_grid<-unique(dft[,c(1,64)]); length(poly_grid$polyID)

# mean cpue --------------------------------------------------------------------

## we have 100 simulated reps for each lat, long, year, season, region and poly
## we summarise the mean cpue for this (need to remove geometry and add it later)
library(plyr)
dftot<-ddply(dft, c("polyID", "Centrelon", "Centrelat", "pseudoyear", "year", 
                    "seas", "yearseas", "region", "psyearID","yearID","regionID"), 
             summarise, mcpue = mean(cpue))

## join to add geometry again
data<-left_join(dftot, poly_grid, by="polyID")

## check
summary(data) ## no NAs
sum(data$mcpue==0) ## no zeros

# save datasets ----------------------------------------------------------------

save.image(file= paste0(out_00,("/00 prepare data.RData")))
rm("df","dft","dftot","lseas","lyear","MyGrid","nyears","xy");ls()
save.image(file="./data/datasets.RData") ## only prepared datasets

# exploratory ------------------------------------------------------------------

## ggpairs
library(GGally)
ggpairs(data[,c(12,6,11,9)])
## save
ggsave(paste0 (out_00,"/00 ggpairs.png"), dpi=300, width=6, height=6)

## histogram continuous vars
library(inspectdf)
data %>% inspect_num %>% show_plot
## save
ggsave(paste0 (out_00,"/00 hist continuous vars.png"), dpi=300, width=6, height=4)

## obs by year
ggplot(data=data, aes(x=as.numeric(as.character(year)), y=mcpue, group=1))+
   geom_col() + xlab("Year") + ylab("Mean CPUE (100 reps)")
## save
ggsave(paste0(out_00, "/00 obs mean cpue by year.png"), dpi=300, width=6, height=4)

## hist by year
qplot(data$mcpue, geom="histogram", fill=data$year, col="black") + 
   xlab("mcpue") + theme_light() + scale_color_manual(values="black")+
   guides(col = FALSE) + labs(fill='Year') 
## save
ggsave(paste0 (out_00, "/00 hist by year.png"), dpi=300, width=6, height=4)

# study area -------------------------------------------------------------------

library(sf)
library(ggspatial)
library(rnaturalearth)
library(maps)
library(ggplot2)

## world map 
world_map_data <- ne_countries(scale = "medium", returnclass = "sf")

## map grid and observed data
ggplot()+
   geom_sf(data = world_map_data) +
   geom_sf(data=data,aes(geometry=geometry), fill="lightblue", alpha=0.5) +
   geom_point(data=data,aes(x=Centrelon, y=Centrelat, col=regionID))+
   theme_light() +
   coord_sf(expand=FALSE, xlim = c(0, 120), ylim = c(-50,50))+
   xlab("longitude") + ylab("latitude")#+   labs(fill='mCPUE')

## save
ggsave(paste0(out_00, "/00 study area poly grids and region pts.png"), dpi=300, width=6, height=6)

## map data frame obs cpue
ggplot() + 
   geom_sf(data = world_map_data) +
   geom_sf(aes(geometry=geometry, fill = mcpue), data = data) +
   scale_fill_viridis_c(option="H") + 
   guides(fill = guide_colourbar(barwidth = 1,
                                 barheight = 10)) +
   theme_light()+
   coord_sf(expand=FALSE, xlim = c(0, 120), ylim = c(-50,50))+
   xlab("longitude") + ylab("latitude") +   labs(fill='mCPUE')

## save
ggsave(paste0(out_00, "/00 study area poly grids and mean cpue pts.png"), dpi=300, width=6, height=6)

#export LD_LIBRARY_PATH=/usr/local/lib
setwd("~/Downloads/cms_analysis")
#password for VPN access: Fratell!n0

#set working directory
#setwd("D/PhD data/Analysis/Chapter3_12Dec17")

####CODE SUMMARY####
##This piece of code takes the raw CMS output to pull out required data corresponding to the specific time steps of interest. This data includes location, exit code, longitude, latitude, depth and distance travelled. The output is saved as a .csv file for each monsoon and season
##The folder organization is as raw_data/*year*/*monsoon* - the subsetted data from the traj_files is saved in this location
##Also creating a list of all the traj_file subsets into a single R data file and saving it to the appropriate raw_data/*year*/*monsoon* location 
####

#install the required packages
#install.packages(c("raster", "sp", "rgdal", "maps", "rgeos", "dplyr", "Hmisc", "ggplot2", "devtools", "digest", "rJava", "geosphere", "stringr", "ncdf4", "sf", "reshape2", "MASS", "spatstat", "maptools", "lattice", "RColorBrewer", "readr", "tidyverse", "tibble", "network", "colorspace", "igraph", "ggraph", "abind", "broom", "viridis"), dependencies=TRUE)

library(raster)
library(sp)
library(rgdal)
library(maps)
library(rgeos)
library(dplyr)
library(Hmisc)
library(ggplot2)
library(devtools)
library(digest)
library(rJava)
library(geosphere)
library(stringr)
library(ncdf4)
library(sf)
library(reshape2)
library(devtools)
library(MASS)
library(spatstat)
library(maptools)
library(lattice)
library(RColorBrewer)
library(readr)
library(ggplot2)
library(tidyverse)
library(network)
library(tibble)
library(colorspace)
library(igraph)
library(ggraph)
library(abind)

gc()

#choosing the required columns from the raw output file of CMS

yr<-c("2009", "2010", "2011")
mon<-c("ne", "sw")

for(g in 1:length(yr)){
for(h in 1:length(mon)){
file.loc<-"~/connectivity-modeling-system/cms-master/expt/cms_output/"
#The following file location is from the blue Seagate
#file.loc<-"G:/backup_CMS_17May19/cms_output/"
#reading the trajectory files based on the year and season information present in the name. g refers to the year and h to the season
traj.all<-list.files(path=paste0(file.loc, yr[g], "/", mon[h]), pattern="*.nc")

#choose the PLD values of interest - this takes the last time step of each day and according to the current PLD range of interest it ranges from 2-50 days with an interval of 2 days in between
#this is for the first time slot of the pld - the first three hours pld.day<-(8*(seq(from=1,to=49, by=2)))+1
pld.day<-8*(seq(from=2,to=50, by=2))

for(i in 1:length(traj.all)){
traj_filename<-paste0(file.loc, yr[g], "/", mon[h], "/", traj.all[i])

#open trajectory file
ncid<-nc_open(traj_filename)

#view the summary of contents in the trajectory file
print(ncid)

#get the values of depth
depth<-ncvar_get(ncid, 'depth')[pld.day,]

#get the values of Longitude
lon<-ncvar_get(ncid, 'lon')[pld.day,]

#get the values of Latitude
lat<-ncvar_get(ncid, 'lat')[pld.day,]

#get the values of exitcode
exit<-ncvar_get(ncid, 'exitcode')

#get the release location id - line number in the release file
loc<-ncvar_get(ncid, 'location')

#get the values of cumulative distance traveled
dis<-ncvar_get(ncid, 'distance')[pld.day,]

#close nestfile
nc_close(ncid)

##writing the data to disk as individual csv files
dat<-rbind(loc, exit, lon, lat, depth, dis)
dat<-t(dat)

colnames(dat)<-c("rel.id", "exit", paste0("lon_", pld.day), paste0("lat_", pld.day), paste0("depth_", pld.day), paste0("distance_", pld.day))

write.csv(dat, paste0("raw_data/", yr[g], "/", mon[h], "/man_op_", substr(traj.all[i], start=1, stop=26), ".csv"), row.names=FALSE)
}
}
}

#save all the man_op files as R objects, so that they are easily accessible
dat.all<-list()

for(g in 1:length(yr)){
for(h in 1:length(mon)){
#list all the trajectory output subset files for a given year and season, corresponding to different release dates
filename<-list.files(path=paste0("raw_data/", yr[g], "/", mon[h]), pattern="man_op_traj_file*")
#save all the files across different release dates in a given year and season as different elements of a list
for(i in 1:length(filename)){
dat.all[[i]]<-read.csv(paste0("raw_data/", yr[g], "/",  mon[h], "/",  filename[i]))
}
#saving dat.all as an R object file - this has all the required PLD subsetted outputs for a given year and season
saveRDS(dat.all, paste0("results/", yr[g], "/", mon[h], "/", yr[g], "_", mon[h], "_man_dat.all.rds"))
}
}

#export LD_LIBRARY_PATH=/usr/local/lib
setwd("~/Downloads/cms_analysis")

#set working directory
#setwd("D/PhD data/Analysis/Chapter3_12Dec17")
  
####CODE SUMMARY####
##This code finds the coordinates of community boundaries across all release date and year replicates for a PLD class per season and uses three methods of plotting these community breaks - the first plots number of breaks per polygon the other two use kernel density estimates of community breaks 
##The output contains -
#1. R data file of a list which saves cut points as a list by monsoon, and lists of PLD classes within each monsoon - this is saved as im_breaks in the location results/infomap_nclus
#2. GeoTiff raster files and jpeg files of distribution of cut points for each monsoon-PLD class combination saved in results/infomap_nclus/
####

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
library(MASS)
library(spatstat)
library(maptools)
library(lattice)
library(RColorBrewer)
library(readr)
library(tidyverse)
library(network)
library(tibble)
library(colorspace)
library(igraph)
library(ggraph)
library(abind)
library(Rcpp)
library(broom)
library(viridis)

gc()

#############################-INPUT DATA-###############################

#setting projection
wgs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#mercator projection
epsg.3395<-"+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#input the polygon file
#I modified the _5Jun18 polygons file to remove gaps that existed between neighbouring polygons to create the _23Oct19 shape file
pol<-readOGR(dsn="dataframes/hycom_small_polygons_final_final_5Jun18", layer="hycom_small_polygons_final_final_10Apr20")

#removing some extra polygons in the Gulf of Khambat
pol2<-pol[-c(62:67),]

pol.large<-readOGR(dsn="dataframes/hycom_polygons_final_final_5Jun18", layer="hycom_polygons_final_final_5Jun18")

#there is an extra polygon in pol.large which needs to be removed and the IDs have to be changed
pol.large<-pol.large[-458, ]

new_IDs<-1:length(pol)
for (i in 1:length(slot(pol.large, "polygons"))){
slot(slot(pol.large, "polygons")[[i]], "ID")<-as.character(new_IDs[i])
}

pol.large$ID<-new_IDs

#input the release points to polygons match file
rel.pol<-read.csv(file="dataframes/rel_small_pol_match.csv", header=TRUE)

#input coastline file - look at readme.txt in all_buf_line_wgs_mod_26Jul19 for steps used to create this coastline file - it says there are 22 features in this file
#I used QGIS to create two files - all_buf_line_wgs_final_26Jul19 and all_buf_line_wgs_22Oct19(where coastline passing through junctions between three-way intersecting polygons are treated differently - look at the QGIS project - all_buf_line_progress which is in the dataframes folder). In all_buf_line_wgs_22Oct19, I have also added line segments joining unconnected nodes in circular island line segments - to me it seems like all the problems with hanging nodes and breaks falling outside polygons have been fixed
#Look at https://gis.stackexchange.com/questions/207945/merging-two-segments-with-overlapping-nodes-in-qgis?rq=1 for solution about removing floating points - snap line segments between selected nodes, merge line segments - copy attributes of the main spatial lines object to the merged object, delete the intervening vertex if necessary
all.buf.line.wgs<-readOGR(dsn="dataframes/all_buf_line_wgs_final_26Jul19", layer="all_buf_line_wgs_final_22Oct19")

#note on 18Feb21 - the command below had worked just fine all the times I ran it before. But now when I run it, it gives me the error - "Error in checkSlotAssignment(object, name, value) : assignment of an object of class “character” is not valid for slot ‘proj4string’ in an object of class “Spatial”; is(value, "CRS") is not TRUE" There is very little documentation of this on the internet, except for the following - http://r-sig-geo.2731867.n2.nabble.com/projection-and-proj4string-td5115767.html. I am unable to understand why the error shows up. I thought it might be because of a possible extra space in the description of the crs, but even when I tried changing it, I received the same error. 
#However, I find that there is a gUnion (works for polygons) analogue for spatial lines called gLineMerge, which can be used. Alternately, I have also merged the features in QGIS and saved the single feature SpatialLines object as all_buf_line_wgs_final_18Feb21. When using this, it seems to detect the same breaks as the combined features operations run using the function aggregate.
#all.buf.line.wgs<-raster::aggregate(all.buf.line.wgs)
all.buf.line.wgs<-gLineMerge(all.buf.line.wgs)

#choose the PLD values of interest - this takes the last time step of each day
#pld.day<-8*(seq(from=2,to=50, by=2))

#I need information on the number of particles released per source polygon for calculating transfer probabilities in the code below - obtaining this from the rel.pol object
rel.pol.count<-as.data.frame(table(rel.pol$rel.pol.id))
colnames(rel.pol.count)<-c("rel.pol.id", "num.part")
rel.pol.count$rel.pol.id<-as.numeric(as.character(rel.pol.count$rel.pol.id))

yrs<-c("2009", "2010", "2011")
mons<-c("ne", "sw")

##################-DETECTING COMMUNITY BOUNDARIES-######################
#input the community detection results that were saved to disk
sum.infomap.all<-read.csv(file="results/infomap_nclus/sum_infomap_all.csv", header=TRUE)
sum.infomap.all$season<-as.character(sum.infomap.all$season)

#creating a vector to save names of plds
pld.sub.names<-c("2_4", "6_12", "14_20", "22_50")

####-REMOVING REDUNDANT CUT POINTS-####

#There are some points within polygons (because of vagaries related to how these polygons were created and how the coastline object runs through them) that do not correspond to cut points but do show up - this gives a false signal of where the cuts are. Creating a dataframe of these points so that they can be removed from the cut points during analysis.
#combining all polygons 
pol.all<-gUnaryUnion(pol)

#intersecting the coastline object with the above
all.int<-gIntersection(all.buf.line.wgs, pol.all, byid=T)

#look at intercept lines one at a time - corresponding to one community at a time
intercept.ends<-coordinates(all.int)
cut.temp.temp<-list()
for(j in 1:length(intercept.ends[[1]])){
start.pt<-intercept.ends[[1]][[j]][1,]
end.pt<-intercept.ends[[1]][[j]][nrow(intercept.ends[[1]][[j]]),]
#For islands which are circular polygons - if all the polygons within the island have the same community membership, there is still a cut point - which is both the start and the end point of the community polygon - removing such points from cut points
if(!identical(start.pt,end.pt)){
cut.temp.temp[[j]]<-as.data.frame(rbind(start.pt, end.pt))
}
}
cut.all<-Reduce(rbind, cut.temp.temp)

####-DETECTING COMMUNITY CUT POINTS-####

#load the rho results
rho.infomap.all<-read.csv(file="results/infomap_nclus/rho_infomap_all.csv", header=TRUE)

#find coordinates of community breaks using community detection results 
for(e in 1:length(yrs)){
cut.mon<-vector("list", length(mons))
for(f in 1:length(mons)){
#select the first season
cut.pld.sub<-list()
for(g in 1:length(pld.sub.names)){
#select the pld class
sum.im.pld<-subset(sum.infomap.all, year==yrs[e] & season==mons[f] & pld.class %in% pld.sub.names[g], select=c(1:534))
rho.im.pld<-subset(rho.infomap.all, year==yrs[e] & season==mons[f] & pld.class %in% pld.sub.names[g], select=c(5:532))

cut.pts<-list()
#go through details of community membership for each row in a PLD class

#creating a dataframe to store the results for each sum.im.pld
res.cut.pol<-sum.im.pld
res.cut.pol[,7:534]<-NA

for(h in 1:nrow(sum.im.pld)){
mem<-t(sum.im.pld[h, 7:ncol(sum.im.pld)])
mem.pol<-as.data.frame(cbind(rel.pol.count$rel.pol.id, mem))
colnames(mem.pol)<-c("pol.id", "infomap.grp")
mem.freq<-as.data.frame(table(mem))
mem.freq$rho<-as.numeric(rho.im.pld[h,1:sum.im.pld$clusters[h]])

#I'll run the code while considering where the membership was greater than a single polygon and also rho is greater than 0.5
mem.freq<-mem.freq[mem.freq$Freq>1,]
mem.freq<-mem.freq[mem.freq$rho>0.5,]

#there are some NA rows generated for some reason, using complete.cases to remove them
mem.freq<-mem.freq[complete.cases(mem.freq),]

#create a list where each element refers to a community
temp<-list()

#take each community for a given simulation and combine the member polygons into one large polygon
for(i in 1:nrow(mem.freq)){
pol.ids<-mem.pol$pol.id[mem.pol$infomap.grp==mem.freq$mem[i]]
temp[[i]]<-gUnaryUnion(pol[pol.ids,])
temp[[i]]@polygons[[1]]@ID<-as.character(i)
}

#create a Spatial polygons object with all the community polygons
temp<-SpatialPolygons(lapply(temp, function(x) slot(x, "polygons")[[1]]))
crs(temp)<-wgs

#code below from - https://gis.stackexchange.com/questions/154689/how-to-find-the-coordinates-of-the-intersection-points-between-two-spatiallines
#intersect community polygons one by one with the coastline Spatial Lines object
intercept<-gIntersection(all.buf.line.wgs, temp, byid=T)

#intercept is a SpatialLines object which has the same number of elements as temp - because intercept is basically the line that runs inside each of the community polygons

cut.temp<-list()
#look at intercept lines one at a time - corresponding to one community at a time
for(i in 1:length(intercept)){
intercept.ends<-coordinates(intercept[i,])
cut.temp.temp<-list()
for(j in 1:length(intercept.ends[[1]])){
start.pt<-intercept.ends[[1]][[j]][1,]
end.pt<-intercept.ends[[1]][[j]][nrow(intercept.ends[[1]][[j]]),]
#For islands which are circular polygons - if all the polygons within the island have the same community membership, there is still a cut point - which is both the start and the end point of the community polygon - removing such points from cut points
if(!identical(start.pt,end.pt)){
cut.temp.temp[[j]]<-as.data.frame(rbind(start.pt, end.pt))
}
}

#EARLIER: If there are two separate polygons for a community, I am taking the northern most and southern most points from that only here - I can choose to consider all of them as well.
#NOW: Since for the circular island coastline - it is hard to separate the gaps from the intercept unless I actually plot it - I am counting each break twice - since one boundary is shared by two polygons
#int<-Reduce(rbind, cut.temp.temp)
#cut.temp[[i]]<-rbind(int[1,], int[length(int),])
#save all the endpoints corresponding to one community, over i such communities from a simulation

cut.raw<-Reduce(rbind, cut.temp.temp)

#removing the extra points within the polygons (which are not cut points from this dataframe) using the anti_join function
if(is.data.frame(cut.raw)){
cut.temp[[i]]<-anti_join(cut.raw, cut.all, by=c("x", "y"))
}
}

#save all end points coming from communities from all simulations of a certain PLD class
cut.pts.pol<-Reduce(rbind, cut.temp)
cut.pts[[h]]<-unique(Reduce(rbind, cut.temp))
print(nrow(cut.pts.pol))
  
#checking if there are duplicated points
#duplicated(cut.pts.pol)
#there are duplicate points in the above - probably because polygons share edges and a cut point falls into both polygons
#I really need to check where these duplicated points are coming from - because very often it is repeated more than once
print(paste0("unique=", nrow(cut.pts.pol[duplicated(cut.pts.pol),])))

#removing duplicates from cut.pts.pol
cut.pts.pol<-unique(cut.pts.pol)

#converting cut.pts.pol into a spatial points object
coordinates(cut.pts.pol)<-~x+y
crs(cut.pts.pol)<-wgs
cut.pol<-over(pol2, cut.pts.pol)
res.cut.pol[h, 7:534]<-ifelse(is.na(cut.pol), 0, 1)   
}

#save dataframe res.cut.pol to disk
write.csv(res.cut.pol, paste0("results/infomap_nclus/infomap_", yrs[e], "_", mons[f], "_", pld.sub.names[g], ".csv"), row.names=FALSE)

#save end points from all PLD classes
cut.pld.sub[[g]]<-Reduce(rbind, cut.pts)
}
cut.mon[[f]]<-cut.pld.sub
#save(cut.pld.sub, file = paste0("results/infomap_nclus/", mons[f], "_im_breaks_1Apr20.RData"))
}
}

#loading all the infomap_ files from the infomap_nclus folder and saving it as a single file
for(i in 1:length(yrs)){
for(j in 1:length(mons)){
for(k in 1:length(pld.sub.names)){
if(i==1 & j==1 & k==1){
cuts.all<-read.csv(file=paste0("results/infomap_nclus/infomap_", yrs[i], "_", mons[j], "_", pld.sub.names[k], ".csv"), header=TRUE)
}
else{
temp<-read.csv(file=paste0("results/infomap_nclus/infomap_", yrs[i], "_", mons[j], "_", pld.sub.names[k], ".csv"), header=TRUE)
cuts.all<-rbind(cuts.all, temp)
}
}
}
}
write.csv(cuts.all, "results/infomap_nclus/infomap_cuts_all.csv", row.names=FALSE)

##################-SAVING POLYGON FILES AS A CSV FILE-################
#using pol2 which has the Gulf of Khambat polygons removed
pol.csv<-fortify(pol2)
str(pol.csv)
pol.csv$id<-as.numeric(pol.csv$id)
pol.csv$id<-pol.csv$id+1

#removing other columns that are not necessary
head(pol.csv)
pol.csv<-pol.csv[,-c(4,5,7)]

colnames(pol.csv)[4]<-"ID"

#saving the dataframe to disk
write.csv(pol.csv, "dataframes/polygons_1Apr20.csv", row.names=FALSE)
###

#export LD_LIBRARY_PATH=/usr/local/lib
#R CMD BATCH filename.r
setwd("~/Downloads/cms_analysis")

#set working directory
#setwd("D/PhD data/Analysis/Chapter3_12Dec17")

####CODE SUMMARY####
##This code reads the R data files saved in raw_data/*year*/*monsoon* one at a time and calculates the transfer probability matrix, converts into a network and runs the Infomap algorithm and saves the community identity output to disk.
##It also pools connectivity matrices for a certain PLD class across release dates and years and writes that to disk as a R data file
##I also calculate connectance values for difference PLD classes across seasons and save it to disk
##the outputs from this code are - 
#1. sum_infomap_*year*_*monsoon*.csv and sum_fgreedy_*year*_*monsoon*.csv (summary of infomap and fastgreedy results) in the location results/infomap_nclus
#2. sum_infomap_all.csv and sum_fgreedy_all.csv which combines the results across years and seasons saved in results/infomap_nclus
#3. man_summary2_*year*_*monsoon*.RData - which saves the connectivity matrices, transfer probability matrices and other miscellaneous output as a R data file in the location results/*year*/*monsoon*
#4. con_mat.RData which saves the connectivity matrices according to PLD class for each monsoon
#5. The last part of the code includes analysis of change in number of clusters with PLD, and statistical tests which gave rise to the PLD classes that were eventually used
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
library(DirectedClustering)
    
gc()

#############################-INPUT DATA-###############################

#setting projection
wgs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#input the polygon file
#pol<-readOGR(dsn="dataframes/hycom_small_polygons_final_final_5Jun18", layer="hycom_small_polygons_final_final_5Jun18")

#I modified the _5Jun18 polygons file to remove gaps that existed between neighbouring polygons
pol<-readOGR(dsn="dataframes/hycom_small_polygons_final_final_5Jun18", layer="hycom_small_polygons_final_final_10Apr20")

pol.large<-readOGR(dsn="dataframes/hycom_polygons_final_final_5Jun18", layer="hycom_polygons_final_final_5Jun18")

#there is an extra polygon in pol.large which needs to be removed and the IDs have to be changed, which is being done in the code below
pol.large<-pol.large[-458, ]

new_IDs<-1:length(pol)
for (i in 1:length(slot(pol.large, "polygons"))){
slot(slot(pol.large, "polygons")[[i]], "ID")<-as.character(new_IDs[i])
}

pol.large$ID<-new_IDs

#input the release points to polygons match file
rel.pol<-read.csv(file="dataframes/rel_small_pol_match.csv", header=TRUE)

#choose the PLD values of interest - this takes the last time step of each day. In the current run, we are using PLDs from 2-50 days
pld.day<-8*(seq(from=2,to=50, by=2))

#I need information on the number of particles released per source polygon for calculating transfer probabilities in the code below - obtaining this from the rel.pol object
rel.pol.count<-as.data.frame(table(rel.pol$rel.pol.id))
colnames(rel.pol.count)<-c("rel.pol.id", "num.part")
rel.pol.count$rel.pol.id<-as.numeric(as.character(rel.pol.count$rel.pol.id))

#from the statistical tests and the distribution of PLDs (which has been commented out), we have chosen four groupings of PLDs - 2-4, 6-12, 14-20, 22-50. Saving these sequences separately as a list 
pld.sub<-list()
pld.sub[[1]]<-c(2,4)
pld.sub[[2]]<-seq(6,12, by=2)
pld.sub[[3]]<-seq(14,20, by=2)
pld.sub[[4]]<-seq(22,50, by=2)

#finding the mean of each of the above PLD ranges for plotting purposes
pld.sub.mean<-sapply(pld.sub, mean)

#saving the different PLD class names for saving outputs
pld.sub.names<-c("2_4", "6_12", "14_20", "22_50")

yrs<-c("2009", "2010", "2011")
mons<-c("ne", "sw")

##################-WRITING OUT CENTRALITY FUNCTIONS-###################
#create a function to calculate in and out strength
#since in and out strength/degree will be plotted on the same graph (as total), and shaded based on one of the components (here it is in), removing self-recruitment from the calculation
in.strength<-function(conn){
in_strength<-colSums(conn)-diag(conn)
}

out.strength<-function(conn){
out_strength<-rowSums(conn)-diag(conn)
}

#create a function to calculate in degree and out degree
in.degree<-function(conn){
conn[conn!=0]<-1
in_degree<-colSums(conn)-diag(conn)
}

out.degree<-function(conn){
conn[conn!=0]<-1
out_degree<-rowSums(conn)-diag(conn)
}

#create a function for calculating clustering coefficient as in the tattoo toolbox
clus.coef<-function(conn){
A<-conn
A[A!=0]<-1
S<-(conn)^(1/3) + (t(conn))^(1/3)
K<-rowSums(A+t(A))
cyc3<-(diag(S %*% S %*% S))/2
K[which(cyc3==0)]<-Inf
cyc3.denom<-(K*(K-1))-(2*(diag(A %*% A)))
C<-cyc3/cyc3.denom
return(C)
}

#create a function to calculate briding centrality (as written by Katell - bridging.m)
bridging<-function(conn){
#converting input connectivity matrix into a binary matrix
BM=conn;
BM[BM!=0]<-1

#calculating total degree
in_deg<-colSums(BM)
out_deg<-rowSums(BM)
Deg<-in_deg+out_deg

#creating an empty vector the size of number of nodes
S<-rep(NA, times=nrow(conn))

for(i in 1:nrow(BM)){
#find indices of nodes to which ith node is connected
v<-which(BM[i,]==1)

#calculating the total degree of connected nodes, inverting it and summing this value across all such nodes
S[i]<-sum(1/Deg[v])
}
Brid<-(1/Deg)/(S)
}

#############-CONNECTIVITY MATRIX AND COMMUNITY DETECTION-#############

#create a list of lists where the outer list is simulation year and the inner lists are connectivity matrices belonging to a certain PLD range that are taken to be replicates of each other across release date.

for(f in 1:length(yrs)){
yr<-yrs[f]

for(g in 1:length(mons)){
mon<-mons[g]

dat.all<-readRDS(paste0("results/", yr, "/", mon, "/", yr, "_", mon, "_man_dat.all.rds"))

##using the final position of each particle for a given pld, running 'over' on these particles with the polygons and finding which particles fall where, convert this into a connectivity matrix

pld.con.mat<-vector("list",12)
pld.rel.od<-vector("list", 12)
pld.tp.all<-vector("list",12)
pld.tp.sim<-vector("list",length(seq(2,50, by=2)))
part.sum<-list()

sum.infomap<-list()
rho.infomap<-list()

bet.centrality<-list()
close.centrality.in<-list()
close.centrality.out<-list()
str.centrality.in<-list()
str.centrality.out<-list()
transv.centrality<-list()
deg.centrality.in<-list()
deg.centrality.out<-list()
brid.centrality<-list()

for(h in 1:length(dat.all)){
#looped over 12 simulations for a monsoonal season
#this loop is to store transfer probabilites across 12 simulations per season and year. In each simulation, 25 PLDs are considered - therefore there are 25 transfer probability matrices in each of the 12 lists in pld.tp.all

#taking the first trajectory file from the list. There are 12 files in the list dat.all in total
dat<-dat.all[[h]]

loc<-dat[,1]
exit<-dat[,2]
lon<-dat[,3:27]
lat<-dat[,28:52]

temp.part<-as.data.frame(matrix(nrow=length(pld.day), ncol=10))
colnames(temp.part)<-c("year", "season", "sim", "pld", "move", "polygon", "sea", "NAs", "land", "out")

for(i in 1:length(pld.day)){
#combine the lon and lat vectors and cbind them
loc.pld<-as.data.frame(cbind(lon[,i], lat[,i], loc, exit))
colnames(loc.pld)<-c("lon", "lat","rel.id","exit")

#finding out how many particles have reached land by exit code -2 and model out of bounds by exitcode -1
loc.pld.na<-loc.pld[is.na(loc.pld$lon),]

#Trajectories could have exit codes -1 and -2 because they exit in the future, but the current locations might be valid - need to include both these statements to ignore trajectories with location NA AND have exit code -2 or -1.

#finding trajectories with coordinates NA and exit code -2 - that end up in land
loc.pld.land<-loc.pld[is.na(loc.pld$lon) & loc.pld$exit==-2,]

#finding trajectories with coordinates NA and exit code -1 - that have gone out of bounds from the model domain
loc.pld.out<-loc.pld[is.na(loc.pld$lon) & loc.pld$exit==-1,]

#keeping rows which have non-NA values
loc.pld.mov<-loc.pld[!is.na(loc.pld$lon),]
#This should be the same as above - loc.pld.mov<-loc.pld[complete.cases(loc.pld$lon),]

#convert it into a spatial points object
coordinates(loc.pld.mov)<-~lon+lat
crs(loc.pld.mov)<-wgs

#run 'over' with loc.pld.mov and the polygon file
set.pol<-over(loc.pld.mov, pol)

#the output has the settlement polygon id as the first column and area of the settlement polygon as the second column, the rownames refer to the release id of the particle
#rownames(set.pol) is identical to loc.pld.mov$rel.id
set.pol<-cbind(rownames(set.pol), set.pol)
colnames(set.pol)<-c("rel.id", "set.pol.id", "set.pol.area")
set.pol$rel.id<-as.numeric(as.character(set.pol$rel.id))
set.pol$set.pol.id<-as.numeric(as.character(set.pol$set.pol.id))

#saving the particles which are still in sea separately. Note that set.pol$rel.id is identical to loc.pld.mov$rel.id - therefore row ids with settlement polygons as NA in set.pol (they haven't reached settlement polygons yet, but are still in sea) will match row ids in loc.pld.mov
loc.pld.sea<-loc.pld.mov[is.na(set.pol$set.pol.id),]

#retaining only those particles which have reached settlement polygons (set.pol$set.pol.id==NA)
set.pol<-set.pol[!is.na(set.pol$set.pol.id),]

#carrying out a left_join using rel.pol and set.pol - keeps all rows from set.pol and and columns from both set.pol and rel.pol. Done to have release and settlement polygon ids side by side.
rel.set<-left_join(set.pol, rel.pol, by="rel.id")

#checking to see that there are no NAs here
which(!complete.cases(rel.set))
colnames(rel.set)

#rel.set is in the format of the con_file output generated by CMS. Using the code I'd written before to include polygons which do not have any settlement

#The spatial polygons dataframe consists of - @data ($ID and $Area), @polygons (list of 534 - $Polygons@ID)
pol.ord<-as.numeric(as.character(pol@data$ID))
npol<-length(pol.ord)

rel<-as.data.frame(cbind(pol.ord, rep(0, times=npol)))
colnames(rel)<-c("rel.pol.id", "dum")

set<-as.data.frame(cbind(pol.ord, rep(0, times=npol)))
colnames(set)<-c("set.pol.id", "dum")

con.mat.in<-left_join(rel, rel.set, by="rel.pol.id")
con.mat.in<-full_join(set, con.mat.in, by="set.pol.id")

#doing a count by running the function table - con.mat.in can have multiple unique combinations of rel and set polygons if multiple connections have occurred between them
con.mat<-as.data.frame(table(con.mat.in$rel.pol.id, con.mat.in$set.pol.id))
colnames(con.mat)<-c("rel.pol.id", "set.pol.id", "count")

con.mat$rel.pol.id<-as.numeric(as.character(con.mat$rel.pol.id))
con.mat$set.pol.id<-as.numeric(as.character(con.mat$set.pol.id))

#checking to see if anything has changed between rel.set and con.mat because of the intervening steps. What we have done with con.mat is added extra combinations of polygons where there was no exchange - adding extra zeros to be exact. Plotting the histogram of non-zero count values between rel.set and con.mat
rel.set.count<-as.data.frame(table(rel.set$rel.pol.id, rel.set$set.pol.id))
colnames(rel.set.count)<-c("rel.pol.id", "set.pol.id", "count")
#test1<-rel.set.count[rel.set.count$count>0,]
#test2<-con.mat[con.mat$count>0,]
#identical(test1$count, test2$count)

#converting the long form to wide form matrix - the release polygons are rows and the settlement polygons are the columns (checked this multiple times)
con.mat.wide<-acast(con.mat, rel.pol.id~set.pol.id, value.var="count")

#h refers to the simulation number (12 simulations per season per year) and i refers to the pld (25 2-day interval PLDs, from 2-50 days)
pld.con.mat[[h]][[i]]<-con.mat.wide

#saving the release id separately, for deriving oceanographic distance. Again - columns: rel.id, exit, lon(columns 3-27), lat(columns 28-52), depth(columns 53-77), distance(columns 78-102)
temp.od<-dat[rel.set$rel.id, c(1, 2+i, 27+i, 52+i, 77+i)]
temp.od$year<-rep(yr, times=nrow(temp.od))
temp.od$mon<-rep(mon, times=nrow(temp.od))
temp.od$sim<-rep(h, times=nrow(temp.od))
temp.od$pld<-rep(pld.day[i]/8, times=nrow(temp.od))

colnames(temp.od)[1:5]<-c("rel.id", "lon", "lat", "depth", "distance")
pld.rel.od[[h]][[i]]<-temp.od

#saving some statistics to disk - the monsoon, the time step within the pld being considered, number of particles on the move, number of particles that have reached a settlement polygon, number of particles still in sea, total number of particles which have exited the simulation, out of the previous - number of particles that have reached land and finally number of particles that have moved out of the model domain. temp.part has 8 rows and there are 4 temp.parts in part.sum
temp.part[i,]<-c(yr, mon, h, pld.day[i]/8, nrow(loc.pld.mov), nrow(set.pol), nrow(loc.pld.sea),nrow(loc.pld.na), nrow(loc.pld.land), nrow(loc.pld.out))

##I need to remove polygons which did not release any particles because, when I divide values by the number of particles released from this - it becomes some number/0 giving -Inf values
pld.con.mat2<-con.mat.wide[-c(62:67), -c(62:67)]

##converting values in pld.con.mat2 into transfer probabilities
#dividing each column in the connectivity matrix with the area of that recipient polygon
pld.tp.mat<-pld.con.mat2
for(j in 1:ncol(pld.tp.mat)){
col.id<-rel.pol.count$rel.pol.id[j]
pld.tp.mat[,j]<-(pld.tp.mat[,j])/(pol@data$Area[pol@data$ID==col.id])
}

#take each row of the above matrix (within the list) and divide by the total number of particles released by the source polygon
for(j in 1:nrow(pld.tp.mat)){
row.id<-rel.pol.count$rel.pol.id[j]
pld.tp.mat[j,]<-(pld.tp.mat[j,])/(rel.pol.count$num.part[j])
}

###NORMALIZING THE CONNECTIVITY MATRIX###
#divide each row by the row sum, so that the sum of all row elements is 1
#for(j in 1:nrow(pld.tp.mat)){
#if(sum(pld.tp.mat[j,])>0){
#pld.tp.mat[j,]<-(pld.tp.mat[j,])/sum(pld.tp.mat[j,])
#}
#}

#write the connectivity matrix to disk - for Katell to cross-check
write.csv(pld.tp.mat, paste0("results/infomap_nclus/con_mat/tp_", yr, "_", mon, "_", "sim", h, "_", pld.day[i]/8, ".csv"), row.names=FALSE)

pld.tp.sim[[i]]<-pld.tp.mat

#clear pld.tp.mat from memory because we will create another matrix with the same name in the loop below
rm(pld.tp.mat)
}

##average values across the pld classes in pld.tp.sim - perform community detection and centrality measure calculations
#create some empty objects where data will be input
temp.infomap<-as.data.frame(matrix(nrow=length(pld.sub), ncol=8+nrow(rel.pol.count)))
colnames(temp.infomap)<-c("year", "season", "sim", "pld.class", "clusters", "clusters2", rel.pol.count$rel.pol.id, "con", "con.prop")

#note that the most number of communities that can be detected are the same as the total number of polygons (i.e. each polygon is a community)
temp.rho<-as.data.frame(matrix(nrow=length(pld.sub), ncol=4+nrow(rel.pol.count)))
colnames(temp.rho)<-c("year", "season", "sim", "pld.class", 1:length(rel.pol.count$rel.pol.id))

#saving centrality measures, the column names are the same as temp.bet, so using the same skeleton
temp.bet<-as.data.frame(matrix(nrow=length(pld.sub), ncol=4+nrow(rel.pol.count)))
colnames(temp.bet)<-c("year", "season", "sim", "pld.class", rel.pol.count$rel.pol.id)
temp.close.in<-temp.bet
temp.close.out<-temp.bet
temp.str.in<-temp.bet
temp.str.out<-temp.bet
temp.transv<-temp.bet
temp.deg.in<-temp.bet
temp.deg.out<-temp.bet
temp.brid<-temp.bet

for(i in 1:length(pld.sub)){
#selecting the PLD sub-classes and converting it into an array of matrices
temp<-abind(pld.tp.sim[pld.sub[[i]]/2], along=3)

#calculating element wise mean across matrices in the array
pld.tp.mat<-apply(temp, c(1,2), mean)

#write the connectivity matrix to disk - for us to cross-check later
write.csv(pld.tp.mat, paste0("results/infomap_nclus/con_mat/tp_", yr, "_", mon, "_", "sim", h, "_", pld.sub.names[i], ".csv"), row.names=FALSE)

#saving pld.tp.mat to the larger list
pld.tp.all[[h]][[i]]<-pld.tp.mat

#calculating connectance values the number and proportion non-zero values in a matrix
con<-length(which(pld.tp.mat>0))
con.prop<-length(which(pld.tp.mat>0))/(length(pld.tp.mat))

#converting the transfer probabilities from wide form to long form for network analysis, also keeping only those rows where the count is greater than 0
pld.tp.long<-melt(pld.tp.mat)
colnames(pld.tp.long)<-c("from", "to", "weight")
pld.tp.long<-pld.tp.long[pld.tp.long$weight>0,]

#Running network analysis with the Infomap algorithm
#creating a directed network with edge weights - note that edge weights are taken as an edge attributes
edges<-pld.tp.long
nodes<-tibble(id=rel.pol.count$rel.pol.id)
routes_igraph<-graph_from_data_frame(d=edges, vertices=nodes, directed=TRUE)

#Infomap algorithm
set.seed(9584)
res.infomap<-cluster_infomap(routes_igraph, e.weights=E(routes_igraph)$weight, modularity=FALSE)

mem<-as.vector(membership(res.infomap))
nclus<-length(unique(mem))
mem.table<-as.data.frame(table(mem))
nclus2<-length(which(mem.table$Freq>1))

temp.infomap[i,]<-c(yr, mon, h, pld.sub.names[i], nclus, nclus2, mem, con, con.prop)

##calculate coherence ratio
#multiplying the conn ectivity matrix with the number of particles released by source polygon to remove that scaling - NEED TO THINK ABOUT THIS A LITTLE BIT
con.scaled<-pld.tp.mat

for(j in 1:nrow(con.scaled)){
con.scaled[j,]<-(con.scaled[j,])*(rel.pol.count$num.part[j])
}

#calculating the total number of particles that were successfully transferred to any other polygon from a given polygon - achieved by calculating row sums for each matrix
con.transfer<-as.data.frame(cbind(rel.pol.count$rel.pol.id, rowSums(con.scaled)))
colnames(con.transfer)<-c("pol.id", "part.success")

com<-as.data.frame(cbind(rel.pol.count$rel.pol.id, mem))
colnames(com)<-c("pol.id", "infomap.grp")
rho<-rep(NA, max(com$infomap.grp))

for(j in 1:max(com$infomap.grp)){
com.pol<-as.character(com$pol.id[com$infomap.grp==j])
#selecting rows and columns based on their names, not row and column numbers
com.mat<-con.scaled[com.pol, com.pol]
com.denom<-con.transfer[com.pol, 2]
rho[j]<-sum(com.mat)/sum(com.denom)
#sometimes the denominator can be 0, if there has been no successful transfer from the polygons in the infomap group - this gives rise to NaN values
}

#saving the rho values as a vector
temp.rho[i,1:(4+length(rho))]<-c(yr, mon, h, pld.sub.names[i], rho)

#building a network where the edge weights are transformed
edges.path<-edges
edges.path$weight<-log10(1/edges$weight)
routes_igraph_path<-graph_from_data_frame(d=edges.path, vertices=nodes, directed=TRUE)

#using the matrix directly instead of using the network made out of it for centrality calculations
routes_igraph_adj<-pld.tp.mat
routes_igraph_path_adj<-log10(1/pld.tp.mat)

#betweenness centrality from the tattoo toolbox
#bet.cen<-betweenness.mod(log10(1/pld.tp.mat))

#number of steps required to reach any other vertex from a given vertex
close.cen.in<-closeness(routes_igraph_path, v=V(routes_igraph_path), weights=E(routes_igraph_path)$weight, mode="in", normalized=TRUE)
temp.close.in[i,]<-c(yr, mon, h, pld.sub.names[i], close.cen.in)

close.cen.out<-closeness(routes_igraph_path, v=V(routes_igraph_path), weights=E(routes_igraph_path)$weight, mode="out", normalized=TRUE)
temp.close.out[i,]<-c(yr, mon, h, pld.sub.names[i], close.cen.out)

#sum of edge weights in and out
str.in<-in.strength(pld.tp.mat)
temp.str.in[i,]<-c(yr, mon, h, pld.sub.names[i], str.in)

str.out<-out.strength(pld.tp.mat)
temp.str.out[i,]<-c(yr, mon, h, pld.sub.names[i], str.out)

#number of edges in and out
deg.in<-in.degree(pld.tp.mat)
temp.deg.in[i,]<-c(yr, mon, h, pld.sub.names[i], deg.in)

deg.out<-out.degree(pld.tp.mat)
temp.deg.out[i,]<-c(yr, mon, h, pld.sub.names[i], deg.out)

#clustering coefficient code as translated from tattoo toolbox in matlab to R by me
FagioloClust<-clus.coef(pld.tp.mat)
temp.transv[i,]<-c(yr, mon, h, pld.sub.names[i], FagioloClust)

#bridging centrality - using the original connectivity matrix to calculate this
brid.cen<-bridging(pld.tp.mat)
temp.brid[i,]<-c(yr, mon, h, pld.sub.names[i], brid.cen)
}
part.sum[[h]]<-temp.part
sum.infomap[[h]]<-temp.infomap
rho.infomap[[h]]<-temp.rho
#bet.centrality[[h]]<-temp.bet
close.centrality.in[[h]]<-temp.close.in
close.centrality.out[[h]]<-temp.close.out
str.centrality.in[[h]]<-temp.str.in
str.centrality.out[[h]]<-temp.str.out
deg.centrality.in[[h]]<-temp.deg.in
deg.centrality.out[[h]]<-temp.deg.out
transv.centrality[[h]]<-temp.transv
brid.centrality[[h]]<-temp.brid
}
#save the connectivity matrices corresponding to each of the PLD classes separately - combining results from all years together into the same list
part.sum<-Reduce(rbind, part.sum)
sum.infomap<-Reduce(rbind, sum.infomap)
rho.infomap<-Reduce(rbind, rho.infomap)
#bet.centrality<-Reduce(rbind, bet.centrality)
close.centrality.in<-Reduce(rbind, close.centrality.in)
close.centrality.out<-Reduce(rbind, close.centrality.out)
str.centrality.in<-Reduce(rbind, str.centrality.in)
str.centrality.out<-Reduce(rbind, str.centrality.out)
deg.centrality.in<-Reduce(rbind, deg.centrality.in)
deg.centrality.out<-Reduce(rbind, deg.centrality.out)
transv.centrality<-Reduce(rbind, transv.centrality)
brid.centrality<-Reduce(rbind, brid.centrality)

write.csv(sum.infomap, paste0("results/infomap_nclus/sum_infomap_", yr, "_", mon, ".csv"), row.names=FALSE)

write.csv(rho.infomap, paste0("results/infomap_nclus/rho_infomap_", yr, "_", mon, ".csv"), row.names=FALSE)

#write.csv(bet.centrality, paste0("results/infomap_nclus/bet_cen_", yr, "_", mon, ".csv"), row.names=FALSE)

write.csv(close.centrality.in, paste0("results/infomap_nclus/close_cen_in_", yr, "_", mon, ".csv"), row.names=FALSE)

write.csv(close.centrality.out, paste0("results/infomap_nclus/close_cen_out_", yr, "_", mon, ".csv"), row.names=FALSE)

write.csv(str.centrality.in, paste0("results/infomap_nclus/str_cen_in_", yr, "_", mon, ".csv"), row.names=FALSE)

write.csv(str.centrality.out, paste0("results/infomap_nclus/str_cen_out_", yr, "_", mon, ".csv"), row.names=FALSE)

write.csv(deg.centrality.in, paste0("results/infomap_nclus/deg_cen_in_", yr, "_", mon, ".csv"), row.names=FALSE)

write.csv(deg.centrality.out, paste0("results/infomap_nclus/deg_cen_out_", yr, "_", mon, ".csv"), row.names=FALSE)

write.csv(transv.centrality, paste0("results/infomap_nclus/trn_cen_", yr, "_", mon, ".csv"), row.names=FALSE)

write.csv(brid.centrality, paste0("results/infomap_nclus/brd_cen_", yr, "_", mon, ".csv"), row.names=FALSE)

#pld.rel.od<-Reduce(rbind, temp.od)
save(part.sum, pld.con.mat, pld.rel.od, pld.tp.all, sum.infomap, rho.infomap, bet.centrality, close.centrality.in, close.centrality.out, str.centrality.in, str.centrality.out, deg.centrality.in, deg.centrality.out, transv.centrality, brid.centrality, file = paste0("results/", yr, "/", mon, "/man_summary_50_",yr, "_", mon, ".RData"))
}
}

##################-PUTTING THE VARIOUS INFOMAP amd RHO OUTPUTS TOGETHER-######################

#loading all the sum_infomap_ files from the infomap_nclus folder and sacing it as a single file

for(i in 1:length(yrs)){
for(j in 1:length(mons)){
if(i==1 & j==1){
sum.infomap.all<-read.csv(file=paste0("results/infomap_nclus/sum_infomap_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
rho.infomap.all<-read.csv(file=paste0("results/infomap_nclus/rho_infomap_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
#bet.all<-read.csv(file=paste0("results/infomap_nclus/bet_cen_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
close.in.all<-read.csv(file=paste0("results/infomap_nclus/close_cen_in_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
close.out.all<-read.csv(file=paste0("results/infomap_nclus/close_cen_out_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
str.in.all<-read.csv(file=paste0("results/infomap_nclus/str_cen_in_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
str.out.all<-read.csv(file=paste0("results/infomap_nclus/str_cen_out_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
deg.in.all<-read.csv(file=paste0("results/infomap_nclus/deg_cen_in_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
deg.out.all<-read.csv(file=paste0("results/infomap_nclus/deg_cen_out_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
trn.all<-read.csv(file=paste0("results/infomap_nclus/trn_cen_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
brd.all<-read.csv(file=paste0("results/infomap_nclus/brd_cen_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
}
else{
temp<-read.csv(file=paste0("results/infomap_nclus/sum_infomap_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
temp.rho<-read.csv(file=paste0("results/infomap_nclus/rho_infomap_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
#temp.bet<-read.csv(file=paste0("results/infomap_nclus/bet_cen_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
temp.close.in<-read.csv(file=paste0("results/infomap_nclus/close_cen_in_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
temp.close.out<-read.csv(file=paste0("results/infomap_nclus/close_cen_out_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
temp.str.in<-read.csv(file=paste0("results/infomap_nclus/str_cen_in_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
temp.str.out<-read.csv(file=paste0("results/infomap_nclus/str_cen_out_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
temp.deg.in<-read.csv(file=paste0("results/infomap_nclus/deg_cen_in_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
temp.deg.out<-read.csv(file=paste0("results/infomap_nclus/deg_cen_out_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
temp.trn<-read.csv(file=paste0("results/infomap_nclus/trn_cen_", yrs[i], "_", mons[j], ".csv"), header=TRUE)
temp.brd<-read.csv(file=paste0("results/infomap_nclus/brd_cen_", yrs[i], "_", mons[j], ".csv"), header=TRUE)

sum.infomap.all<-rbind(sum.infomap.all, temp)
rho.infomap.all<-rbind(rho.infomap.all, temp.rho)
#bet.all<-rbind(bet.all, temp.bet)
close.in.all<-rbind(close.in.all, temp.close.in)
close.out.all<-rbind(close.out.all, temp.close.out)
str.in.all<-rbind(str.in.all, temp.str.in)
str.out.all<-rbind(str.out.all, temp.str.out)
deg.in.all<-rbind(deg.in.all, temp.deg.in)
deg.out.all<-rbind(deg.out.all, temp.deg.out)
trn.all<-rbind(trn.all, temp.trn)
brd.all<-rbind(brd.all, temp.brd)
}
}
}

#save the output to disk
write.csv(sum.infomap.all, "results/infomap_nclus/sum_infomap_all.csv", row.names=FALSE)
write.csv(rho.infomap.all, "results/infomap_nclus/rho_infomap_all.csv", row.names=FALSE)
#write.csv(bet.all, "results/infomap_nclus/bet_raw.csv", row.names=FALSE)
write.csv(close.in.all, "results/infomap_nclus/close_in_raw.csv", row.names=FALSE)
write.csv(close.out.all, "results/infomap_nclus/close_out_raw.csv", row.names=FALSE)
write.csv(str.in.all, "results/infomap_nclus/str_in_raw.csv", row.names=FALSE)
write.csv(str.out.all, "results/infomap_nclus/str_out_raw.csv", row.names=FALSE)
write.csv(deg.in.all, "results/infomap_nclus/deg_in_raw.csv", row.names=FALSE)
write.csv(deg.out.all, "results/infomap_nclus/deg_out_raw.csv", row.names=FALSE)
write.csv(trn.all, "results/infomap_nclus/trn_raw.csv", row.names=FALSE)
write.csv(brd.all, "results/infomap_nclus/brd_raw.csv", row.names=FALSE)

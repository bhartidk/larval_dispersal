setwd("D:/PhD data/Analysis/Chapter3_12Dec17")

library(lubridate)

#input coordinates file
rel.loc<-read.csv(file="dataframes/release_wgs_final_22May18.csv")
nrow(rel.loc)

colnames(rel.loc)<-c("Longitude", "Latitude")

rel.loc$Polygon=rep(1, times=nrow(rel.loc))
rel.loc$Depth=rep(1, times=nrow(rel.loc))
rel.loc$Number=rep(1, times=nrow(rel.loc))

#first date of release
year<-2010
month<-7
day<-1
hour<-0
mins<-0
secs<-0

#creating columns for the time of release information
rel.loc$Year<-rep(year, times=nrow(rel.loc))
rel.loc$Month<-rep(month, times=nrow(rel.loc))
rel.loc$Day<-rep(day, times=nrow(rel.loc))
rel.loc$Second<-rep(secs, times=nrow(rel.loc))

#looking at the resultant dataframe
head(rel.loc)

#reorganizing the columns
rel.loc<-rel.loc[,c(3, 1, 2, 4:ncol(rel.loc))]
head(rel.loc)

rel.date<-ymd_hms(paste(year, month, day, hour, mins, secs, sep="-"))

ndays<-10
time.steps<-ndays*8
rel.loc.fin<-rel.loc
rel.loc.int<-rel.loc

for(i in 1:(time.steps-1)){
day1<-rel.date
hour(day1)<-3*i
date1<-c(year(day1), month(day1), day(day1), (hour(day1)*3600))
rel.loc.int[,6:ncol(rel.loc.int)]<-t(replicate(nrow(rel.loc.int), date1))
rel.loc.fin<-rbind(rel.loc.fin, rel.loc.int)
}

nrow(rel.loc.fin)
head(rel.loc.fin)
#getting the end date required for nest_1.nml file and the timeMax specification required for the runconf.list file

#outputing the dataframe as a tab delimited text file
#write.table(rel.loc.fin, file=paste0("dataframes/releaseFiles/releaseFile_", day, "_", month, "_", year), sep=" ", row.names=FALSE, col.names=FALSE)

write.table(rel.loc.fin, file=paste0("dataframes/releaseFiles/releaseFile_",day, "_", month, "_", year), sep=" ", row.names=FALSE, col.names=FALSE)

#last release date
last.rel<-rel.date +  hours((time.steps-1)*3)
last.rel

#last run date - this would be 120 days + last.rel date
#specifying timeMax
timeMax.days<-120

#timeMax in seconds - timeMax.days*24 hours*60 minutes*60 seconds
timeMax<-timeMax.days*24*60*60
timeMax

# 
last.run<-last.rel + days(timeMax.days)
last.run

#checking if the number of rows matches the number of rows there should be there considering 2136 locations * 8 time steps per day * 10 days - it matches!
nrow(rel.loc.fin)

#checking to see if the time-steps have been correctly assigned - looking at the unique columns corresponding to date in the release file
rel.disp<-unique(rel.loc.fin[,6:ncol(rel.loc.fin)])
rel.disp


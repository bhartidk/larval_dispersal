%loading the different identifiers for files
year={'2009','2010','2011'};
season={'ne','sw'};
sim={'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'};
PLDclass={'2_4','6_12','14_20','22_50'};

%creating numeric values corresponding to the identifiers above
year_names=2009:2011;
season_names=1:2;
sim_names=1:12;
PLDclass_names=1:4;

%creating an empty matrix where betweenness centrality measures will be
%stored
bet_all=[];

%loading files by year, season, sim and pld class
for nyear=1:3
for nseason=1:2
for nsim=1:12
for npldclass=1:4

%put together the file name that needs to be loaded
file=strcat('tp_',year{nyear},'_',season{nseason},'_sim', sim(nsim), '_', PLDclass(npldclass),'.csv');

%loading the file
con=readmatrix(file{1});

%removing the first row which is the header
con=con(2:529,:);

%converting the weights into path lengths by using newproba, because
%betweenness involves the calculation of shortest paths
con=newproba(con);

%calculating betweenness centrality
bet=betweenness_w(con);
bet=bet';

%saving betweenness centrality with some of the other identifiers
id=[nyear, nseason, nsim, npldclass];
bet_all=[bet_all; id, bet];
end
end
end
end

%save results as a csv file
writematrix(bet_all,'bet_mat.csv') 
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

%creating an empty matrix where subgraph summary will be
%stored
sub_sum=[];

%creating an empty matrix where closeness centrality will be stored
sub_close=[];

%create a vector with nodenames - 61-68 would be missing
nodes_g=[1:61 68:534];
nodes_g=string(nodes_g');

%nodes_g=cellstr(num2str(nodes_g'));
%nodes_g=convertCharsToStrings(nodes_g);

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

%converting con into a graph object
con_g=digraph(con, nodes_g, 'omitselfloops');

%find the number of connected graphs
[bin,binsize]=conncomp(con_g, 'Type', 'weak');

%save the total number of singleton nodes
sub_single=sum(binsize==1);
sub_ids=find(binsize>1);
sub_multi=length(sub_ids);
sub_msize=mean(binsize(sub_ids));

%saving subgraph summary
id=[nyear, nseason, nsim, npldclass];
sub_sum=[sub_sum; id, sub_single, sub_multi, sub_msize];

end
end
end
end

%save results as a csv file
writematrix(sub_sum,'subgraph_weak_sum.csv'); 

%creating a function for calculating bridging coefficient
function Brid = bridging(conn)
%converting input connectivity matrix into a binary matrix
BM=conn~=0;

%calculating total degree
[in_deg,out_deg,Deg]=degrees_dir(BM);

%creating another empty matrix
S=[];

%looping over the number of nodes

for k=1:size(BM,1)

%find indices of nodes have non-zero connectivity values
v=find(BM(k,:)==1);

%calculating the total degree of connected nodes, inverting it and summing this value across all suc nodes
S=[S,sum((1./Deg(v)))];

end
Brid=(1./Deg)./(S);
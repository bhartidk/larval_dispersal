tic

clear all

datein = [datenum(2009,1,9,0,0,0):datenum(0,0,0,3,0,0):datenum(2009,6,1,0,0,0)];
dateout = datestr(datein,'yyyymmddHHMMSS');
s1 = 'nest_1_';
s2 = '.nc';
file = strcat(s1, dateout, s2);

%%% path of the folder which has the input data and output data folders
%path = strcat(pwd, '/');
%path = '/home/bharti/Downloads/cms_input/2010_ne/';
%path = '/home/bharti/connectivity-modeling-system/cms-master/expt/expt_extra/download_raw/';
path = '/home/bharti/Seagate/backup_CMS_19Dec18/input/hycom_2009_raw/';

%%% name of the input data folder
data = '';

%%% name of the output data folder
%new_data = 'test/';
%new_data = '/home/nlp/Documents/';
new_data = '/home/bharti/connectivity-modeling-system/cms-master/expt/expt_2010/nests/raw/';

nc_names = strcat(path, data, file);

for inc = 1:(size(nc_names,1)-1)

%%% open nestfile
t1 = nc_names(inc,:);
t2 = nc_names(inc+1,:);
nc1 = netcdf.open(t1,'NOWRITE');
nc2 = netcdf.open(t2,'NOWRITE');

%%%% extract variable from netcdf file
%%%% variables requested
%%% surf_el = free surface elevation at time t   => eta(i,j)
%%% surf_el = free surface elevation at time t+1   => eta_f(i,j)
%%% water_u = horizontal velocity in eastward => u(i,j,k) 
%%% water_v = horizontal velocity in northward => v(i,j,k)
%%% lon = longitude of mesh grid
%%% lat= latitude of mesh grid

varid_lon=netcdf.inqVarID(nc1, 'lon');
longi=netcdf.getVar(nc1,varid_lon);

varid_lat=netcdf.inqVarID(nc1, 'lat');
lati=netcdf.getVar(nc1,varid_lat);

varid_z=netcdf.inqVarID(nc1, 'depth');

%%% z is positive downward in HYCOM outputs but we need it positive upward for w calculation
z=-netcdf.getVar(nc1,varid_z);   

varid_t=netcdf.inqVarID(nc1, 'time');
t=netcdf.getVar(nc1,varid_t);

varid_sal=netcdf.inqVarID(nc1, 'salinity');
scale_factor = ncreadatt(t1,'salinity','scale_factor');
add_offset = ncreadatt(t1,'salinity','add_offset');
sal=netcdf.getVar(nc1,varid_sal);
index=find(sal==-3e+04);
sal(index)=0;
sal=double(sal)*scale_factor+add_offset;

varid_temp=netcdf.inqVarID(nc1, 'water_temp');
scale_factor = ncreadatt(t1,'water_temp','scale_factor');
add_offset = ncreadatt(t1,'water_temp','add_offset');
temp=netcdf.getVar(nc1,varid_temp);
index=find(temp==-3e+04);
temp(index)=0;
temp=double(temp)*scale_factor+add_offset;

varid_u=netcdf.inqVarID(nc1, 'water_u');
scale_factor = ncreadatt(t1,'water_u','scale_factor');
add_offset = ncreadatt(t1,'water_u','add_offset');
u=netcdf.getVar(nc1,varid_u);
index=find(u==-3e+04);
u(index)=0;
u=double(u)*scale_factor+add_offset;

varid_v=netcdf.inqVarID(nc1, 'water_v');
scale_factor = ncreadatt(t1,'water_v','scale_factor');
add_offset = ncreadatt(t1,'water_v','add_offset');
v=netcdf.getVar(nc1,varid_v);
index=find(v==-3e+04);
v(index)=0;
v=double(v)*scale_factor+add_offset;

varid_eta=netcdf.inqVarID(nc1, 'surf_el');
scale_factor = ncreadatt(t1,'surf_el','scale_factor');
add_offset = ncreadatt(t1,'surf_el','add_offset');
eta=netcdf.getVar(nc1,varid_eta);
eta_index=find(eta==-3e+04);
eta(eta_index)=0;
eta=double(eta)*scale_factor+add_offset;

varid_eta_f=netcdf.inqVarID(nc2, 'surf_el');
scale_factor = ncreadatt(t2,'surf_el','scale_factor');
add_offset = ncreadatt(t2,'surf_el','add_offset');
eta_f=netcdf.getVar(nc2,varid_eta);
etaf_index=find(eta_f==-3e+04);
eta_f(etaf_index)=0;
eta_f=double(eta_f)*scale_factor+add_offset;

[imax,jmax,kmax]=size(u);

%% time step dt in seconds
dt=3*3600;%%%
R=6371000;
%% horizontal grid steps dx and dy calculated after Mercator projection of lon and lat
for i=2:imax-1
for j=2:jmax-1
dx(i,j)=(longi(i+1)-longi(i-1))/2/360*2*pi*R*cos(lati(j)/180*pi);
dy(i,j)=(lati(j+1)-lati(j-1))/2/360*2*pi*R;
end
end

%%%% caution: k levels are also increasing downward (k=0 is free surface and k=kmax is bottom)
%%%% so defining dz = z(k-1) - z(k), ensures dz positive upward => ensures w positive upward 
dz=z(1:kmax-1)-z(2:kmax);

%%% compute free surface vertical velocity
for i=2:imax-1
for j=2:jmax-1
%%%% compute vertical velocity at the free surface  => ws
w(i,j,1)=(eta_f(i,j)-eta(i,j))/dt+u(i,j,1)*(eta(i+1,j)-eta(i-1,j))/(2*dx(i,j))+v(i,j,1)*(eta(i,j+1)-eta(i,j-1))/(2*dy(i,j));
end
end

%%% compute vertical velocity along depth
for k=2:kmax
for i=2:imax-1
for j=2:jmax-1
%%%% compute vertical velocity at the free surface  => ws
w(i,j,k)=w(i,j,k-1)+dz(k-1)*((u(i+1,j,k)-u(i-1,j,k))/(2*dx(i,j))+(v(i,j+1,k)-v(i,j-1,k))/(2*dy(i,j)));
end
end
end

%% create a new nc file with the modified variables 

%%% the following line suppresses scaling that matlab applies on a
%%% matrix if it has extreme values http://kb.mit.edu/confluence/pages/viewpage.action?pageId=3907432
format long g;

%%% scaling all the variables extracted above 
sf_all=1e-3;  %% scale factor for u, v, salinity, temperature, surface elevation
sf_w=1e-6;  %% scale factor for w

sal_int=int32(sal/sf_all);
temp_int=int32(temp/sf_all);
u_int=int16(u/sf_all);
v_int=int16(v/sf_all);
eta_int=int16(eta/sf_all);
w_int=int32(round(w/sf_w));

%%% changing the sign on depth levels so that it reads positive downwards to
%%% make it CMS compatible
z=-z;
%%% changing the sign on w-velocity so that it reads positive downwards to
%%% make it CMS compatible
w_int=-w_int;

%%% define the fill value for u, v and surface elevation to the one recognized by CMS
fill_val=int16(-3e+04);
%%% define fill value for salinity and temperature separately which is in int32 format 
fill_val_st=int32(-3e+04);
%%% define fillvalue for w will be recognized as -30 after applying the scale
%%% factor 1e-6 on it
fill_val_w=int32(-3e+7);

%%% apply fill value to land mask points
sal_int(index)=fill_val_st;
temp_int(index)=fill_val_st;
u_int(index)=fill_val;
v_int(index)=fill_val;
eta_int(eta_index)=fill_val;

%%% find the indices where w is 0 - land mask to apply fill value
w_index=find(w==0);
w_int(w_index)=fill_val_w;

%%% create a new nc file with the same id as the first file read in
%%% nc_new ='nest_1_test.nc';
nc_new = strcat(path, new_data, file(inc,:));

%%% NC_NOCLOBBER to prevent overwriting of any file with the same name
nc3 = netcdf.create(nc_new,'NC_NOCLOBBER');

%%% define new dimensions
dimidtime = netcdf.defDim(nc3,'time',length(t));
dimidlat = netcdf.defDim(nc3,'lat',length(2:jmax-1));
dimidlon = netcdf.defDim(nc3,'lon',length(2:imax-1));

%%% subset by depth
depth_levels = 27;

dimiddepth = netcdf.defDim(nc3,'depth', depth_levels);

%%% define variables
varid_eta = netcdf.defVar(nc3,'surf_el','short',[dimidlon dimidlat dimidtime]);
varid_t = netcdf.defVar(nc3,'time','double',[dimidtime]);
varid_lat = netcdf.defVar(nc3,'lat','double',[dimidlat]);
varid_lon = netcdf.defVar(nc3,'lon','double',[dimidlon]);
varid_sal = netcdf.defVar(nc3,'salinity','int',[dimidlon dimidlat dimiddepth dimidtime]);
varid_z = netcdf.defVar(nc3,'depth','double',[dimiddepth]);
varid_temp = netcdf.defVar(nc3,'water_temp','int',[dimidlon dimidlat dimiddepth dimidtime]);
varid_u = netcdf.defVar(nc3,'water_u','short',[dimidlon dimidlat dimiddepth dimidtime]);
varid_v = netcdf.defVar(nc3,'water_v','short',[dimidlon dimidlat dimiddepth dimidtime]);

%%% create fill value attribute
netcdf.putAtt(nc3,varid_sal,'_FillValue',fill_val_st);
netcdf.putAtt(nc3,varid_sal,'scale_factor',sf_all);

netcdf.putAtt(nc3,varid_temp,'_FillValue',fill_val_st);
netcdf.putAtt(nc3,varid_temp,'scale_factor',sf_all);

netcdf.putAtt(nc3,varid_u,'_FillValue',fill_val);
netcdf.putAtt(nc3,varid_u,'scale_factor',sf_all);

netcdf.putAtt(nc3,varid_v,'_FillValue',fill_val);
netcdf.putAtt(nc3,varid_v,'scale_factor',sf_all);

netcdf.putAtt(nc3,varid_eta,'_FillValue',fill_val);
netcdf.putAtt(nc3,varid_eta,'scale_factor',sf_all);

%%% close netcdf file definitions
netcdf.endDef(nc3);

%%% copy variable data into the netcdf file
netcdf.putVar(nc3, varid_eta, eta_int(2:imax-1,2:jmax-1));
netcdf.putVar(nc3, varid_t, t);
netcdf.putVar(nc3, varid_lat, lati(2:jmax-1));
netcdf.putVar(nc3, varid_lon, longi(2:imax-1));
netcdf.putVar(nc3, varid_sal, sal_int(2:imax-1,2:jmax-1,1:depth_levels,:));
netcdf.putVar(nc3, varid_z, z(1:depth_levels));
netcdf.putVar(nc3, varid_temp, temp_int(2:imax-1,2:jmax-1,1:depth_levels,:));
netcdf.putVar(nc3, varid_u, u_int(2:imax-1,2:jmax-1,1:depth_levels,:));
netcdf.putVar(nc3, varid_v, v_int(2:imax-1,2:jmax-1,1:depth_levels,:));

%% create a new nc file for the extra w velocity variable
%%% saving the w velocity nc file separately
nc_w = strcat(new_data, s1, dateout(inc,:), 'w', s2);
nc4 = netcdf.create(nc_w,'NC_NOCLOBBER');

%%% define new dimensions
dimidtime = netcdf.defDim(nc4,'time',1);
dimidlat = netcdf.defDim(nc4,'lat',length(2:jmax-1));
dimidlon = netcdf.defDim(nc4,'lon',length(2:imax-1));
dimiddepth = netcdf.defDim(nc4,'depth', depth_levels);

%%% define variables
varid_t_w = netcdf.defVar(nc4,'time','double',[dimidtime]);
varid_lat_w = netcdf.defVar(nc4,'lat','double',[dimidlat]);
varid_lon_w = netcdf.defVar(nc4,'lon','double',[dimidlon]);
varid_z_w = netcdf.defVar(nc4,'depth','double',[dimiddepth]);
varid_w = netcdf.defVar(nc4,'water_w','int',[dimidlon dimidlat dimiddepth dimidtime]);

%%% create fill value attribute
netcdf.putAtt(nc4,varid_w,'_FillValue',fill_val_w);
netcdf.putAtt(nc4, varid_w, 'scale_factor', sf_w); 

%%% close netcdf file definitions
netcdf.endDef(nc4);

%%% copy variable data into the netcdf file
netcdf.putVar(nc4, varid_t_w, t);
netcdf.putVar(nc4, varid_lat_w, lati(2:jmax-1));
netcdf.putVar(nc4, varid_lon_w, longi(2:imax-1));
netcdf.putVar(nc4, varid_z_w, z(1:depth_levels));
netcdf.putVar(nc4, varid_w, w_int(2:imax-1,2:jmax-1,1:depth_levels,:));

%close netcdf files
netcdf.close(nc1);
netcdf.close(nc2);
netcdf.close(nc3);
netcdf.close(nc4);
       
end

toc





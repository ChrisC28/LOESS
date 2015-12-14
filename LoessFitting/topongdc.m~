% TOPONGDC:   Extract heights from topo_ngdc datafile, for given lats&lons.
% INPUT:
%    lat - matrix of lats
%    lon - matrix of longs, 0-360deg
%    fname - optional - name of terrainbase file
%
% OUTPUT:
%    deps - in metres, +ve upwards (so actually heights, not depths).
%
% Jeff Dunn 26/5/97
% 22/5/98  Sped up by reading in topo file chunks
%
% USAGE: [deps] = topongdc(lats,lons,fname)

function [deps] = topongdc(lats,lons,fname)

%ncmex('setopts',0);


% Changed from 6.2 to 8.2 on 8/1/2003
if nargin<3 | isempty(fname)
   fname = '/home/netcdf-data/topo_ngdc_8.2.nc';
end



deps = repmat(NaN,size(lats));
idx = 1:prod(size(lats));

length(idx)
ierrcnt = 0;
rcode = 0;
rcode = -1;
%while rcode<0
%   
    %[fid,rcode] = ncmex('ncopen',fname,'nowrite');
%   if rcode<0
%      ierrcnt = ierrcnt+1;
%      if ierrcnt>50
%	 error('Aborting in TOPONGDC - flaky connection to bathy file');
%      else
%	 pause(1);
%      end
%   end
%end

%flon = ncread(fname,'lon');
%flat = ncread(fname,'lat');

load(fname);
flon = topo_lon;
flat = topo_lat;


jj = find(lats<flat(1) | lats>flat(end));
indices_in_domain = setxor(1:length(lats),jj);

if ~isempty(jj)
  disp([num2str(length(jj),8) ' lats out of range']);
  idx(jj) = [];
end

for i=1:length(jj) 
    lats(jj(i))
end 

% 17/3/05 Without properly tracking down the problem, avoid an indexing
% overflow just by avoiding lon=360.
jj = find(lons>=360);
if ~isempty(jj)
   lons(jj) = 359.99;
end

% NOTE: varget uses 'C' style indexing from 0


max_lat = max(lats(indices_in_domain));
min_lat = min(lats(indices_in_domain));
max_lon = max(lons(indices_in_domain));
min_lon = min(lons(indices_in_domain));

acc = 0.5;


index_max_lat = find(flat <= round(max_lat/acc)*acc + acc); 
index_max_lat = index_max_lat(end);

index_min_lat = find(flat >= round(min_lat/acc)*acc - acc); 
index_min_lat = index_min_lat(1);

index_max_lon = find(flon <= round(max_lon/acc)*acc + acc); 
index_max_lon = index_max_lon(end);

index_min_lon = find(flon >= round(min_lon/acc)*acc - acc); 
index_min_lon = index_min_lon(1);


%netcdf_file_id = netcdf.open(fname,'NC_NOWRITE');
%height_var_id  = netcdf.inqVarID(netcdf_file_id,'height');
%height_var_id  = netcdf.inqVarID(netcdf_file_id,'height');
%height = netcdf.getVar( netcdf_file_id, height_var_id, [index_min_lon, index_min_lat],  ... 
%		       [index_max_lon-index_min_lon,index_max_lat-index_min_lat]);

height = topo_depth(index_min_lon:index_max_lon-1,index_min_lat:index_max_lat-1);

deps=interp2(flon(index_min_lon:index_max_lon-1),flat(index_min_lat:index_max_lat-1), height', lons(indices_in_domain), lats(indices_in_domain),'cubic' ); 

%netcdf.close(netcdf_file_id);


%{
lorec = floor(lons(idx).*30);

rlt0 = flat(1);
nlat = length(flat);
nlon = length(flon);
rad = pi/180;
alpha = rad/30; 
beta= tan(rad*(45.+rlt0/2.));
size(lats(idx))

larec = floor(-.5 + (log(tan(rad*(45+lats(idx)./2))/beta))/alpha);
% Reading one netcdf element at a time is slow, so if more than a few points
% required, build a mesh of chunks and read in and extract from any chunks 
% containing required points.

chsz = 500;
netcdf_file_id = netcdf.open(fname,'NC_NOWRITE');
height_var_id  = netcdf.inqVarID(netcdf_file_id,'height');

%dch = ncmex('varget',fid,'height',[lach(ii) loch(ii)],[chs1 chs2]);

if length(idx)>100
  [lach,loch] = meshgrid(min(larec):chsz:max(larec),min(lorec):chsz:max(lorec));
  for ii=1:prod(size(lach))
    jj = find( larec>=lach(ii) & larec<(lach(ii)+chsz) ...
 	     & lorec>=loch(ii) & lorec<(loch(ii)+chsz) );
    if ~isempty(jj)
       if (lach(ii)-1+chsz) > nlat
	  chs1 = nlt-lach(ii);
       else
	  chs1 = chsz;
       end
       if (loch(ii)-1+chsz) > nlon
	  chs2 = nlon-loch(ii);
       else
	  chs2 = chsz;
       end
       %dch = ncmex('varget',fid,'height',[lach(ii) loch(ii)],[chs1 chs2]);
      
       dch   = netcdf.getVar(netcdf_file_id, height_var_id, [lach(ii)+1,loch(ii)+1], [chs1, chs2]);

       latmp = larec(jj)-lach(ii)+1;
       lotmp = lorec(jj)-loch(ii)+1;
       chidx = (chs2*(latmp-1))+lotmp;
       deps(idx(jj)) = dch(chidx);
    end
  end  
else
  for jj = 1:length(idx)
    %deps(idx(jj)) = ncmex('varget1',fid,'height',[larec(jj) lorec(jj)]);
    %larec(jj)
    %lorec(jj) 
    deps(idx(jj)) = netcdf.getVar(netcdf_file_id, height_var_id,[larec(jj) lorec(jj)]);
  end  
end
%}


%ilat_south=find(lats<-70);
%load /usr/home/cchlod/LoessFitting/req_deps.mat
%req_deps(find(req_deps==0))=NaN;
%cpt=0;f
%for ii=1:length(ilat_south);
%   if (lons(ilat_south(ii))<0) llon=lons(ilat_south(ii))+360; else llon=lons(ilat_south(ii)); end
%   dep_c=interp2(rgx,rgy,req_deps, llon, lats(ilat_south(ii))); 
%   if ~isnan(dep_c)
%      deps(ilat_south(ii))=dep_c;
%      cpt=cpt+1;
%   end
%end
%disp(['Recover ' num2str(cpt) ' lats']);


% ------------ End of topongdc.m -------------------

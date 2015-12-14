% GET_BATH  Get ocean depth at given locations, using a range of datasets.
%
% INPUTS
%  x,y    Locations at which depth is required
%  dset   vector of datasets in order
%         1=Terrainbase  2=NGDCv8.2  3=AusBath15   4=AGSO98  5=AGSO2002
%         Default:  [2 5]
%  rpt    [Optional]  1 => save loaded AusBath15 to global variable, to make
%         subsequent calls faster. Default 0.
% OUTPUTS
%  deps   depths (m), -ve downwards, NaN where no value.
%
% SUPERCEDES  get_bath15, agso_bath_xy
%
% Jeff Dunn   CSIRO Marine Research   8/1/2003
%
% SEE ALSO   get_bath_agso.m   (to get AGSO 2002 at full resolution)
%
% USAGE: deps = get_bath(x,y,dset,rpt)
%
%

function deps = get_bath(x,y)
workspace;

% Mods: 9/5/03 Pre-check for presence of datasets, and resort to others
%       if some are missing.

% Adapted by JB SALLEE (2006 to work at legos(Toulouse)
% simplified to only use NGDCv8.2

%ncquiet;

% What sort of computer are we running on?
% Assuming Unix
plat = 0;
repDataperso = num2str('/net/argos/data/peps/cchlod/')
fnm = num2str([repDataperso num2str('/Bathy/ETOPO1_360_wrap_SH.mat')]);
%fnm = num2str('topo_ngdc_8.2.nc');

aa = exist(fnm,'file');

if aa==0
   disp('GET_BATH: Cannot find any of the bathy files - no depths returned.')
   return
end

deps = repmat(nan,size(x));

ii = 1:prod(size(x));
       % NGDC
       deps(ii) = topongdc(y(ii),x(ii),fnm); 
end
   
%---------------------------------------------------------------------------

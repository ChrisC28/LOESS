% REQ_CHECK   return a list of the given grid points which DO need to
%              be mapped, (the rest are too shallow).
% INPUT
%     lon,lat   Coordinates of points to check (any dimension)
%     dep       Mapping depth (m +ve downwards) to check against. May be either
%               scalar or same size as lon
%     opt       grid-spacing=required if any of grid cell corners required 
%               [default 0=required if grid point (or cell centre) required]
%     excl_poly [np 2] polygon defining region to exclude. Default []
%
% Jeff Dunn 6/5/96
%   Mod:   14/5/98  Generalise to arbitrary req_grid, and also keep output in
%                   order.
%          28/3/00  Convert to depths rather than standard depths
%          26/3/02  Allow for new .25 degree grid req_deps.mat
%          11/6/02  Add 'opt', so that instead of checking grid point, check
%                   all 4 corners of grid cell.
%           
% $Id: req_check.m,v 1.10 2002/06/11 04:22:59 dun216 Exp dun216 $
%
% USAGE: req = req_check(lon,lat,dep,opt)

function req = req_check(lon,lat,dep,opt,excl_poly)

persistent REQ_CHECK_alertflag

% Argument 5 used to activate exclusion of points in the North Atlantic.
% Now made more general by allowing a polygon to define region to exclude.

if nargin<5 || isempty(excl_poly) || all(excl_poly(:)==0)
   excl_poly = [];
end

if nargin<4 | isempty(opt)
   opt = 0;
end

load req_deps;

if any(lat>max(rgy(:)))
   if isempty(REQ_CHECK_alertflag)
      disp('REQ_CHECK: latitude north of limit of req_deps fields');
      REQ_CHECK_alertflag = 1;
   end
end

% Initialise so that do not map points outside of 'req_deps'
req = zeros(size(lon));

jj = find(lon(:)<=max(rgx) & lon(:)>=min(rgx) & lat(:)<max(rgy) ...
       & lat(:)>min(rgy));

dep = dep + .1;      % Increase a fraction so that 0m doesn't map on land, etc
if max(size(dep))>1
   dep = dep(:);
   djj = jj;
else
   djj = 1;
end

if opt~=0
   ox = [-opt opt opt -opt]/2;
   oy = [opt opt -opt -opt]/2;
   for ii=1:4
      tmp = req_chk(lon(jj)+ox(ii),lat(jj)+oy(ii),dep(djj),req_deps,rginc,rgx,rgy);
      if size(req,2) ~= size(tmp,2)
	 tmp = tmp';
      end
      req(jj) = req(jj) | tmp;
   end
else
   req(jj) = req_chk(lon(jj),lat(jj),dep(djj),req_deps,rginc,rgx,rgy);
end

req = find(req>0);

if ~isempty(excl_poly)
   ii = inpolygon(lon(req),lat(req),excl_poly(:,1),excl_poly(:,2));
   req(ii) = [];
end


return

%---------------------------------------------------------------------------
function req = req_chk(lon,lat,dep,req_deps,rginc,rgx,rgy)

if any(lon(:)<0 | lon(:)>=360)
   ii = find(lon<0);
   lon(ii) = lon(ii)+360;
   ii = find(lon>=360);
   lon(ii) = lon(ii)-360;
end

if all(rem(lon*2,1)==0) & all(rem(lat*2,1)==0) & rginc==.5

  % Old "req_deps" grid was 1/2 degree, so could do efficient direct calc 
  % (2 cols for every degree longitude, + 2 x lat degrees (from 1 grid below
  % min latitude).
  lcol = length(rgy);
  ii = ((lon-min(rgx))*2*lcol) + ((lat-min(rgy)+.5)*2);
  req = req_deps(ii)>=dep;

elseif all(rem(lon*4,1)==0) & all(rem(lat*4,1)==0) & rginc==.25

  % This "req_deps" grid is 1/4 degree, so can do efficient direct calc
  % (4 cols for every degree longitude, + 4 x lat degrees (from 1 grid below
  % min latitude).
  lcol = length(rgy);
  ii = ((lon-min(rgx))*4*lcol) + ((lat-min(rgy)+.25)*4);
  req = req_deps(ii)>=dep;

else

  [xx,yy] = meshgrid(rgx,rgy);
  rg = interp2(xx,yy,req_deps,lon,lat,'*linear');
  req = (rg(:)>=dep);
  
end

%----------------------------------------------------------------------

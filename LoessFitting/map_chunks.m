% MAP_CHUNKS  - apply Bogus Bathy loess mapping by sub-chunks, w/wo BAR, with
%               data pre-weighting (usually downweighting values estimated
%               from adjacent levels using local z gradient.)
% INPUT:  
%  lo, la    Longitude and latitude vectors of entire datasource region
%  doy       Day of year vector (or days-since-1900 if fitting fns 3 or 4)
%  dd        Data vector
%  dw        adjacent-level weight (presquared; 0=min 1=max downweight)
%  bdep      Sea bottom depth at data locations - MUST BE non -ve
%  xG,yG     Map locations (optional NaNs where mapping not required) 
%  gdep      Sea bottom depth at grid locations
%  config    Structure containing (at least) the fields:
%  - depm       Depth of mapping, in m
%  - q1,q2      Loess mapping parameters
%  - rmin,rmax  More loess mapping parameters
%  - fns        [1 n]  n=1:annual  n=2:semian  n=3:T  n=4:T^2  {def [0 0 0 0]}
%  - bar        1 = use BAR system
%  - xw         XWeight - zonal stretch coefficient for loess radius
%  - bthy       0 = do NOT use (TAR); otherwise amount of TAR (eg .7 or 1)
%  - m2coast    1=> map all non-inland locations, irrespective of depth
%               2=> map everywhere, irrespective of land/ocean 
%  - bogus      0=off  1=Spatial  2=Spatial&Temporal -ve=bogus-to-zero
%  - smw        1=return sum of weights  [def 1]
%  - lowrap     1=offset longitudes so that do not need to handle Greenwich wrap
%  - var_rmin   1=use rmin field input argument  [def 0]
%  - var_xw     1=use xw field input argument  [def 0]
%      and may contain:
%  - noextrap   1=do not map where purely extrapolating
%               2=retry mapping with doubled q1 if map fails due to 
%               extrapolation test. Only expect to use for gapfilling by
%               remapping!
%               3,4=do not map if fails more fancy extrapolation tests
%  xm,ym     Range beyond grid that we seek data (km, lat corrected). Empty
%            arguments triggers auto-adaptive chunk extension.
%  rmin      rmin field if variable rmin
%  xw        xw field is variable xw
%
% OUTPUT: mn: Mean
%         rq: distance (km) from xi,yi to the most distant (nq'th) data point'
%         nq: number of points used in each mapping (or <=0 error code)
%   Optional:
%         an: Annual harmonics    (returned empty if not required)
%         sa: Semi-annual harmonics  (returned empty if not required)
%         swgt: sum-of-weights
%         tco: extra fitted function (Time)
%         t2co: extra fitted function (Time^2)
%
% USAGE: [mn,rq,nq,an,sa,swgt,tco,t2co] = map_chunks(lo,la,doy,dd,dw,bdep,...
%        xG,yG,gdep,config,xm,ym,rmin,xw)

function [mn,rq,nq,an,sa,swgt,tco,t2co] = map_chunks(lo,la,doy,dd,dw,bdep,...
				xG,yG,gdep,config,xm,ym,rmin,xw)

% Branched from:   
% $Id: map_chunks.m,v 1.2 2007/05/14 01:10:14 dun216 Exp dun216 $
% Copyright:  Jeff Dunn   CSIRO Marine Research
%
% MODS: 30/3/05  Collect control args into "config", and set up treatment of
%          longitude wrap.

global SHUTUP_MAP_CHUNKS

% Chunk size, in degrees.
xinc = 5; yinc = 5;
%minext = xinc/4;
minext = xinc/3;   % 13/4/06 Changed from xinc/4 prior to CARS2005 2nd pass
disp('hello')
if isempty(config.rmax)
   config.rmax = 10e10;
end
if ~isfield(config,'var_rmin') || ~config.var_rmin
   rmin = config.rmin;
end
if ~isfield(config,'var_xw') || ~config.var_xw
   xw = config.xw;
end
if isempty(config.bthy)
   config.bthy = 1;
end

if isempty(config.fns)
   config.fns = [0 0 0 0];
elseif length(config.fns)==1 & config.fns==2
   config.fns = [1 1 0 0];
elseif length(config.fns)<4      
   config.fns(4) = 0;
end

if isempty(config.bar)
   config.bar = 0;
end
if isempty(config.smw)
   config.smw = 1;
end
if isempty(config.lowrap)
   config.lowrap = 0;
end

if nargin<11
  xm = xinc;
  ym = yinc;
elseif ~isempty(xm)
   % Convert margins from km to degrees
   ym = ym/111;  xm = xm/111;
end
   
if config.bar
   global  P_Gate P_line P_P
   global  lstxw pp_xw_init
end

dd = dd(:);

if isempty(dw)
   dw = zeros(size(dd));
end

if length(xw(:))==1
   xw = repmat(xw,size(xG));
end

if length(rmin(:))==1
   rmin = repmat(rmin,size(xG));
end


start_time = cputime;

if any(isnan(bdep))
   ii = find(isnan(bdep));
   bdep(ii) = get_bath(lo(ii),la(ii));
end
if any(bdep<0)
   ii = find(bdep<0);
   bdep(ii) = 0;
end

tmp = repmat(nan,size(xG));
mn = tmp;
rq = tmp;
nq = tmp;

[an,sa,tco,t2co,swgt] = deal([]);

if config.fns(1);  an = tmp;    end
if config.fns(2);  sa = tmp;    end
if config.fns(3);  tco = tmp;   end
if config.fns(4);  t2co = tmp;  end
if config.smw;   swgt = tmp;    end


% Resolve grid interval so req_check can test grid cell corners
if min(size(yG))<=1
   ginc = 0;
else
   ginc = abs(yG(2)-yG(1));
   if ~ginc
      ginc = abs(xG(2)-xG(1));
   end
end

ii = find(~isnan(xG));
if config.m2coast==2
   req = ii;
elseif config.m2coast==1 || config.depm<0
   req = ii(req_check(xG(ii),yG(ii),0,ginc));
else
   req = ii(req_check(xG(ii),yG(ii),config.depm,ginc));
end
   
if isempty(gdep)
   gdep = zeros(size(req(:)));
else
   gdep = gdep(req(:));
end
xG = xG(req);
yG = yG(req);
rmin = rmin(req);
xw = xw(req);

if config.bar
   % Determine in which polygons the data points and grid points lie.
   % The polygons are stored in a [3 n] cell array called 'poly'. It is loaded 
   % from polys_N.mat (where N=depth level).
   % cells(1,n) scalar Id for n-th poly
   % cells(2,n) vector of lons of corners for n-th poly
   % cells(3,n) vector of lats of corners for n-th poly
   % 
   % Note that polygon with ID 1 is a catchall default, and is not in Poly.

   % Use the deepest completed depth at or above the mapping depth.
   % max([dep 0]) means that if dep<0, which indicates isopycnal mapping,
   % then will use the sea surface BAR
   disp('Do something stupid')
   load('/home/dunn/eez/bar_data/completed');
   ii = find(completed<=max([config.depm 0]));
   pdep = completed(max(ii));

   load(['/home/dunn/eez/bar_data/polys_' num2str(pdep)]);   
   
   pold = ones(size(lo));
   poli = ones(size(xG));
	 
   for ii = 1:size(regionpoly,2)
      isin = find(inpolygon(lo,la,regionpoly{2,ii},regionpoly{3,ii}));
      if ~isempty(isin)
	 pold(isin) = repmat(regionpoly{1,ii},size(isin));
      end
      isin = find(inpolygon(xG,yG,regionpoly{2,ii},regionpoly{3,ii}));
      if ~isempty(isin)
	 poli(isin) = repmat(regionpoly{1,ii},size(isin));
      end
   end

   mxpo = max([pold(:); poli(:)]);

   % Load global poly-relation arrays from file
   ppload(pdep,mxpo);

   if config.lowrap
      % Save copy of BAR so that can switch between 0 and 180E centred longitude
      P_line_save = P_line;
      P_Gate_save = P_Gate;
      P_P_save = P_P;
   end
end


if config.lowrap
   % Save copy of longitudes so that can switch between 0 and 180E centred
   lo_save = lo;
end
wrapchng = 0;
wrapon = 0;

% Make mesh of chunks
[xch,ych] = meshgrid(min(xG(:)):xinc:max(xG(:)),min(yG(:)):yinc:max(yG(:)));
[m n] = size(xch);

maxey = config.rmax/111;

% Process each chunk
for jj = 1:m*n
  GG = find(xG>=xch(jj) & xG<(xch(jj)+xinc) & yG>=ych(jj) & yG<(ych(jj)+yinc));
  if ~isempty(GG)
    lcor = latcor(mean(yG(GG)));
    
    if ~isempty(xm)
       ymin = min(yG(GG))-ym;
       ymax = max(yG(GG))+ym;
       xmin = min(xG(GG))-xm/lcor;
       xmax = max(xG(GG))+xm/lcor;

       if ~config.lowrap
	  % then do not do following tests
       elseif wrapon>=0 & xmax>360
	  % If was 0 or 180 centred but now need to be 360 centred...
	  wrapon = -1; wrapchng = 1;
	  lo = lo_save;
	  kk = find(lo<180);
	  lo(kk) = lo(kk)+360;
       elseif wrapon<=0 & xmin<0
	  % If was 360 or 180 centred but now need to be 0 centred...
	  wrapon = 1; wrapchng = 1;
	  lo = lo_save;
	  kk = find(lo>180);
	  lo(kk) = lo(kk)-360;
       elseif (wrapon==1 && xmax>180) || (wrapon==-1 && xmin<180)
	  % If was 360 or 0 centred but now need to be 180 centred...
	  wrapon = 0; wrapchng = 1;
	  lo = lo_save;
       end

       II = find(lo>=xmin & lo<=xmax & la>=ymin & la<=ymax);

    else
       % 6/5/05  Coded this adaptive chunk margins algorithm
       cbox = [min(xG(GG)) max(xG(GG)) min(yG(GG)) max(yG(GG))];
    
       maxex = config.rmax/(lcor*111);
       obox = cbox + [-maxex maxex -maxey maxey];

       if ~config.lowrap
	  % then do not do following tests
       elseif wrapon>=0 & obox(2)>360
	  % If was 0 or 180 centred but now need to be 360 centred...
	  wrapon = -1; wrapchng = 1;
	  lo = lo_save;
	  kk = find(lo<180);
	  lo(kk) = lo(kk)+360;
       elseif wrapon<=0 & obox(1)<0
	  % If was 360 or 180 centred but now need to be 0 centred...
	  wrapon = 1; wrapchng = 1;
	  lo = lo_save;
	  kk = find(lo>180);
	  lo(kk) = lo(kk)-360;
       elseif (wrapon==1 && obox(2)>180) || (wrapon==-1 && obox(1)<180)
	  % If was 360 or 0 centred but now needs to be 180 centred...
	  wrapon = 0; wrapchng = 1;
	  lo = lo_save;
       end

       II = find(lo>=obox(1) & lo<=obox(2) & la>=obox(3) & la<=obox(4));
       if length(II)>(config.q1*5)
	  Xx = (lo(II)-cbox(1))*lcor;
	  Xx = Xx.*Xx;
	  Yy = la(II)-cbox(3);
	  r = sqrt(Xx + Yy.*Yy);
	  rmx(1) = binsort(r,config.q1);
	  Yy = la(II)-cbox(4);
	  Yy = Yy.*Yy;
	  r = sqrt(Xx + Yy);
	  rmx(2) = binsort(r,config.q1);
	  Xx = (lo(II)-cbox(2))*lcor;
	  Xx = Xx.*Xx;
	  r = sqrt(Xx + Yy);
	  rmx(3) = binsort(r,config.q1);
	  Yy = la(II)-cbox(3);
	  r = sqrt(Xx + Yy.*Yy);
	  rmx(4) = binsort(r,config.q1);
	  
	  % Expand the selected region by 10%, as cutting it too fine 
	  % sometimes leaves block structure in maps
	  rext = max(rmx)*1.1;       
	  if rext<minext
	     rext = minext;
	  end
	  if (rext/lcor)<maxex | rext<maxey
	     tmpx = min([maxex rext/lcor]);
	     tmpy = min([maxey rext]);	  
	     obox = cbox + [-tmpx tmpx -tmpy tmpy];
	     II = find(lo>=obox(1) & lo<=obox(2) & la>=obox(3) & la<=obox(4));
	  end
       end
    end
    
    basin_box=1 ;
    rmoutlier=1;
    if ~isempty(II)
     if basin_box==1;
      load polybasin.mat
      if xch(jj)<0 xch_c=xch(jj)+360; else xch_c=xch(jj); end
      lo_c=lo(II);
      lo_c(find(lo_c<0))=lo_c(find(lo_c<0))+360;
      irm=[];
      if inpolygon( xch_c, ych(jj), xA1, yA)
         irm1=inpolygon(lo_c,la(II), xI, yI);						   
         irm2=inpolygon(lo_c,la(II), xP, yP);
         irm=unique([find(irm1==1) find(irm2==1)]);
      elseif inpolygon( xch_c, ych(jj), xA2, yA)
         irm1=inpolygon(lo_c,la(II), xI, yI);						   
         irm2=inpolygon(lo_c,la(II), xP, yP);
         irm=unique([find(irm1==1) find(irm2==1)]);
      elseif inpolygon( xch_c, ych(jj), xI, yI)
         irm1=inpolygon(lo_c,la(II), xA1, yA);						   
         irm2=inpolygon(lo_c,la(II), xA2, yA);						   
         irm3=inpolygon(lo_c,la(II), xP, yP);
         irm=unique([find(irm1==1) find(irm2==1) find(irm3==1)]);
      elseif inpolygon( xch_c, ych(jj), xP, yP)
         irm1=inpolygon(lo_c,la(II), xA1, yA);						   
         irm2=inpolygon(lo_c,la(II), xA2, yA);						   
         irm3=inpolygon(lo_c,la(II), xI, yI);
         irm=unique([find(irm1==1) find(irm2==1) find(irm3==1)]);
      end
      if ~isempty(irm)
         
         II(irm)=[];
      end
     end
    end				   
    if ~isempty(II)
     if rmoutlier
       mmmm=nanmean(dd(II));
       ssss=nanstd(dd(II));
       irm=find(dd(II)>mmmm+3*ssss | dd(II)<mmmm-3*ssss);
       if ~isempty(irm)
         II(irm)=[];
       end
     end
    end
    Gro = req(GG);
    
    % config carries the following through to loessb and loessbbar:
    %  q1, q2, rmax, fns, bthy, bogus, smw
    if config.bar
       
       % If just swapped longitude centre, need to change all longs in BAR
       % definitions and have BAR distances recalculated.
       if wrapchng==1
	  wrapchng = 0;

	  % Change lstxw and clear the init status so that BAR paths are 
	  % recalculated
	  lstxw = -1;
	  pp_xw_init = zeros(size(pp_xw_init));
	  
	  P_line = P_line_save;
	  P_Gate = P_Gate_save;
	  P_P = P_P_save;
       
	  if wrapon==0
	     % Don't need to modify the retrieved original BAR definitions'
		  
	  elseif wrapon==-1
	     kk = find(P_Gate(1,:)<180);
	     P_Gate(1,kk) = P_Gate(1,kk)+360;

	     kk = find(P_line(1,1,:)<180);
	     P_line(1,1,kk) = P_line(1,1,kk)+360;
	     kk = find(P_line(1,2,:)<180);
	     P_line(1,2,kk) = P_line(1,2,kk)+360;

	     for pp = 1:size(P_P,2)
		for mm = 1:3
		   if ~isempty(P_P{6,pp}{mm})
		      kk = find(P_P{6,pp}{mm}(1,:)<180);
		      P_P{6,pp}{mm}(1,kk) = P_P{6,pp}{mm}(1,kk)+360;
		   end
		end
	     end       
	  elseif wrapon==1
	     kk = find(P_Gate(1,:)>180);
	     P_Gate(1,kk) = P_Gate(1,kk)-360;
	     
	     kk = find(P_line(1,1,:)>180);
	     P_line(1,1,kk) = P_line(1,1,kk)-360;
	     kk = find(P_line(1,2,:)>180);
	     P_line(1,2,kk) = P_line(1,2,kk)-360;

	     for pp = 1:size(P_P,2)
		for mm = 1:3
		   if ~isempty(P_P{6,pp}{mm})
		      kk = find(P_P{6,pp}{mm}(1,:)>180);
		      P_P{6,pp}{mm}(1,kk) = P_P{6,pp}{mm}(1,kk)-360;
		   end
		end
	     end       
	  end
       end
       
       vo = loessbbar(lo(II),la(II),doy(II),bdep(II),dd(II),dw(II),...
		      xG(GG),yG(GG),gdep(GG),pold(II),poli(GG),rmin(GG),...
		      xw(GG),config);


       % Rarely used tacky complication 26 June 07.
       % If filling mid-ocean gaps by remapping (using existing map values
       % as datapoints) then want to bridge gaps rather than give up 
       % because actually extrapolating from both sides. Succeed in a few
       % more cases if retry with q1 doubled (so larger r, hence more chance
       % of see points on other side of gap.)
       if isfield(config,'noextrap') & config.noextrap==2 & any(vo.nq==-4)
	  tmpcfg = config;
	  tmpcfg.q1 = 2*config.q1;
	  i4 = find(vo.nq==-4);
	  G4 = GG(i4);
	  vo2 = loessbbar(lo(II),la(II),doy(II),bdep(II),dd(II),dw(II),...
			  xG(G4),yG(G4),gdep(G4),pold(II),poli(G4),rmin(G4),...
			  xw(G4),tmpcfg);
	  fnms = fieldnames(vo2);
	  for ifn = 1:length(fnms)
	     eval(['vo.' fnms{ifn} '(i4) = vo2.' fnms{ifn} ';']);
	  end
	  
       end

	  else
       vo = loessb(lo(II),la(II),doy(II),bdep(II),dd(II),dw(II),...
		   xG(GG),yG(GG),gdep(GG),rmin(GG),xw(GG),config);
    end
    mn(Gro) = vo.mn;
    rq(Gro) = vo.rq;
    nq(Gro) = vo.nq;
    if config.fns(1)
       an(Gro) = vo.an;
    end
    if config.fns(2)
       sa(Gro) = vo.sa;
    end
    if config.fns(3)
       tco(Gro) = vo.tco;
    end
    if config.fns(4)
       t2co(Gro) = vo.t2co;
    end
    if config.smw
       swgt(Gro) = vo.swgt;
    end
  end
end

if isempty(SHUTUP_MAP_CHUNKS)
   disp(['MAP_CHUNKS cpu time = ' num2str(cputime-start_time,10) ' s']);
end

% ------------ End of map_chunks.m ---------------------

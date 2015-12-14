% LOESS_MAP:  loess map user-supplied data, via map_chunks
%
% PATHS required:
%        /home/toolbox/local/csirolib/
%        /home/dunn/matlab/
%        /home/dunn/matlab/saga/
%        /home/dunn/eez/general/
%        /home/eez_data/software/matlab/
%
% INPUT:
%  xg,yg  Map locations, pref plaid grid (otherwise do NOT use interpfill
%         or screen or collect stats)
%         NOTE: grid should not cross 0E - ie map Atlantic as 300:360 and
%         then 0:20, not as 300:20 nor -60:20.) 
%  xd,yd  Data locations
%  dd     Data values (can be complex)
%  tims   Unix time of data points (days since 1/1/1900)
%  dw     data weights (0-1). Use [] if not required.
%           dw=0   => ignore point
%           dw=.5  => 1/2 weight
%           dw=1   => full weight - fully use data point
%  fns    [1 n]  n=1:annual  n=2:semian  n=3:T  n=4:T^2  {def [0 0 0 0]}
%  bar    1=use BAR (BAR depth level from value of 'depm') [def 0]
%  scrn   0=just do mapping  [default]
%         1=calc stats and do residuals screening   
%         2=just calc stats but not screen
%  fopt   Vector of option codes 
%     1=rmin  Loess min data radius. 0=> use CARS pre-computed field
%         of variable rmin  [def 300]
%     2=rmax  Loess max data radius   [def 1500]
%     3=xw   Loess zonal stretch (>1 => zonal stretch)
%        0 =>  use CARS pre-computed field of variable xw  
%     5=TAR  1=use bathy-weighting (TAR) system  [default 1]
%     11=m2coast  0=only map where ocean is >=depm deep [default 0] 
%                 1=map any ocean locations
%                 2=map anywhere, even if not ocean
%     12=scrthr - multiple of RMSR => screening threshold  [def 3]   
%        Only used if scrn=1
%     13=bogus   [def 1]  0=none  1=spatial
%           2-4=spatial plus increasing levels of temporal bogusing
%           -ve=use a bogus of 0, weight: -1=1, -2=sumw/100 -3=sumw/30
%               -4=temporal bogus to 0
%           (effectively a background field of 0 where data is scarce)
%     14=Q  Base number of data points to use **This will be scaled for
%        the number of fit functions; eg if fns = [1 1] then Q=240
%        will translate to Q=400 actually used. [default 240]
%     15=sumwgts  1=compute sum-of-weights for each map point
%     16=smth_dep 1=presmooth bathymetry before applying TAR, to 
%        reduce topo-induced structure, esp in hi-res coastal maps 
%     18=interpfill  Use locally-weighted interpolating of surrounding
%        mapped values to fill gaps where data distribution inadequate
%        for loess interpolation (only with plaid grid) [def 1]
%        A drawback is that this does NOT use BAR or TAR.
%     19=rmsmin  Min allowed RMSR to multiply by scrthr for a screening
%         limit.
%     23=auto_ext  1=calculate minimum required chunksizes rather than 
%         fixed sizes, for faster mapping, especially where lots
%         of data. May sometimes result in a little block structure.
%     26=noextrap  1= do NOT map if purely extrapolating (set nq=-4) 
%         2= remap with double q1 if fails simple extrapolation test 
%         3,4= do NOT map if fails more subtle extrapolation tests (nq=-5,-6)
%          [default 1]
%     31-34=fn_damp_ranges  - ask Jeff
%     40=comment - description of property being mapped etc, for 
%         details string
%     41=output file name [default - 'lmap_0' if depm=0]
%        0 = do NOT produce output file! 
%     42=depm  [default 0]  Nominal mapping depth, used to determine 
%         BAR level, t-function-dampoff [where used], and default
%         output file name. 
%     44=alternate rmin filename (if provided will be used instead of 
%        normal XW and RMIN fields.)
%     45=wgtest  1=set nq=-3 and do not map if sum-of-wgts too low [def 1]
%
%  fval   cell array of values assigned to the parameters specified by fopt
%
% OUTPUT:     map variables written to mat-file       
%    config - structure containing mapping configuration info
%    details - string describing mapping
%    xg,yg  - if NOT a plaid grid, otherwise vectors Xg,Yg
%    mn     - mean field [see also nq, below]
%    nq     - number of points used for each mapping. **Values <=0 flag
%              pathological
%    rq     - data radius (km) used for each mapping
%  Optionally...
%    an     - complex annual cycle (fns(1)=1)
%    sa     - complex semi-annual cycle (fns(2)=1)
%    tco    - linear time trend (fns(3)=1)
%    t2co   - quadratic of time trend (fns(4)=1)
%    swgt   - sum of weights
%    sd     - locally-weighted standard deviation of datapoints used 
%    rmsmr  - RMS of residuals of data wrt mean field
%    rmsr   - RMS of residuals of data wrt seasonal field (fns 1 or 2 used)
%    rejidx - index to data points which exceed scrthr rmsr
%
%   IF OUTPUT ARGS
%    mn, nq,rq  are returned
%
% Note that the coefficients of the spatial quadratic fit are not saved.
%
% AUTHOR: Jeff Dunn  CSIRO (CMAR) Feb 2007
%   (devolved from lodrive.m)
%
%   Longitude wrap is handled.
%
% EXAMPLE: 
% > [xg,yg] = meshgrid(100:.5:180,-50:.5:0);
%     say we want to down-weight older data, and have data from 1950-2000 
%     arrange weights 1950=0 to 2000=1
% > dw = (year-1950)./50
%     want annual and semi-annual harmonics
% > fns = [1 1]
%     set rmin=300km, zonal stretch=1.5, strong temporal bogussing, leave
%     gaps where loess fails rather than fill, add comment to descriptive
%     text in output file, output to "blah1.mat" :
% > fopt = [1 3 13 18 40 41];
% > fval = {300,1.5,4,0,'my blah data, weights favour recent data','blah1'};
%
% > loess_map(xg,yg,xd,yd,dd,tims,dw,fns,1,0,fopt,fval)
% 
% USAGE: [mn,nq,rq] = loess_map(xg,yg,xd,yd,dd,tims,dw,fns,bar,scrn,fopt,fval);

function [mn,nq,rq] = loess_map(xg,yg,xd,yd,dd,tims,dw,fns,bar,scrn,fopt,fval)

global lstxw pp_xw_init ppname


% Control variables, modifiable by fopt - see lodrive_opts.m
depm = 0;
rmin = 300;
rmax = 1500;
interpfill = 1;
TAR = 1;
xw = 1;
m2coast = 0;
var_rmin = 0;
rmin_fnam = [];
var_xw = 0;
Q = 240;
scrthr = 3.0;
bogus=1;
smw = 1;
smth_dep = 0;
rmsmin = 0;
auto_ext = 1;
noextrap = 1;
damp{4} = [];
rQCn = 1;
cmt = [];
fnm = [];
wgtest = 1;

if ~isempty(dw)
   % Convert from 1=full weight -> 1=max down-weighting, as required by
   % loess routines. Then pre-square, also required.
   dw = (1-dw).*(1-dw);
   dw=dw(:);
end
if nargin < 10 | isempty(scrn);  scrn = 0; end
if nargin < 9 | isempty(bar);  bar = 0; end

if length(fns)==1 & fns==2
   % Convert from old to new style "fns" 
   fns = [1 1 0 0];
elseif length(fns)<4      
   fns(4) = 0;
end


if nargin>=11 & ~isempty(fopt) & ~isempty(fval)
   for ii = 1:length(fopt)
      switch fopt(ii)
	case 1
	  rmin = fval{ii};
	case 2
	  rmax = fval{ii};
	case 3
	  if fval{ii}>0
	     xw = fval{ii};
	     var_xw = 0;
	  else
	     xw = 1;
	     var_xw = 1;
	  end	   
	case 5
	  TAR = fval{ii};
	case 11
	  m2coast = fval{ii};
	case 12
	  scrthr = fval{ii};
	case 13
	  bogus = fval{ii};
	case 14
	  Q = fval{ii};
	case 15
	  smw = fval{ii};
	case 16
	  smth_dep = fval{ii};
	case 18
	  interpfill = fval{ii};
	case 19
	  rmsmin = fval{ii};
	case 23
	  auto_ext = fval{ii};
	case 24
	  rQCn = fval{ii};
	case 26
	  noextrap = fval{ii};
	case {31,32,33,34}
	  damp{fopt(ii)-30} = fval{ii};
	case 40
	  cmt = fval{ii};
	case 41
	  fnm = fval{ii};
	case 42
	  depm = fval{ii};
	case 44
	  rmin_fnam = fval{ii};
	case 45
	  wgtest = fval{ii};
	otherwise
	  disp(['*** Do not understand FOPT = ' num2str(fopt(ii))]);
      end
   end
end

if rmin==0
   var_rmin = 1;
else
   var_rmin = 0;
end
ym = rmax;
xm = floor(ym*xw);



if var_xw | var_rmin
   if ~isempty(rmin_fnam)
      load(rmin_fnam)
   else
      load /home/dunn/eez/map_data/loess_xw
   end
   II = find(~isnan(xg) & xg>xfmin & xg<xfmax & yg>yfmin & yg<yfmax);
   if var_xw
      xw = repmat(xwdef,size(xg)); 
      xw(II) = interp2(xfield,yfield,xw_field,xg(II),yg(II),'*linear');
   end
   if var_rmin
      rmin = repmat(rmindef,size(xg)); 
      rmin(II) = interp2(xfield,yfield,rmin_field,xg(II),yg(II),'*linear');
   end
   clear xfield yfield rmin_field xw_field
end

if bar
   lstxw = -1;
   pp_xw_init = zeros(1,11);
   ppname = ['tmp_pp_' num2str(depm) '_'];
end

if ~isempty(damp)
   % Activate constrained-function loess, with constraint ramped over
   % a user-defined depth range.

   constr = 0;
   ncof = [2 2 1 1];
   % Assume 6 spatial coeffs - ie quadratic plus x*y terms.
   %expow = 8;
   expow = 6;
   ful = 10.^expow;
   si = repmat(ful,1,6);
   for ii = 1:4
      if fns(ii)
	 jj = length(si) + (1:ncof(ii));
	 if ~isempty(damp{ii})
	    if depm<=damp{ii}(1)
	       % leave default 'full-fit' coefficient
	       dexp = expow;
	    elseif depm>=damp{ii}(2)
	       % total damping coefficient
	       disp(['LODRIVE: At or below damp-off depth for function ' ...
		     num2str(ii)]);
	       jj = [];
	       fns(ii) = 0;
	       dexp = 0;
	    else
	       dexp = 2 - 6*(depm-damp{ii}(1))./diff(damp{ii});
	       constr = 1;
	    end
	    si(jj) = 10.^dexp;
	 else	       
	    si(jj) = ful;
	 end
      elseif ~isempty(damp{ii})
	 disp(['LOESS_MAP Warning: Damping specified for non-fitted ' ...
                  'function ' num2str(ii)]);
      end
   end
   if ~constr
      % If above or below the damping range for all functions, then no 
      % constraints required.
      si = [];
   end
else
   si = [];
end


% Set required number of points on basis of functions being fitted. 28/9/98
nfns = sum(fns);
Q = Q*(1 + nfns/3);
q1 = round(Q);
q2 = round(Q/5);
disp(['Adjusting for fit-fns, req number of points: ' num2str([q1 q2])]);



cfns = '';
if all(fns==0); cfns = 'Spatial only';     end
if fns(1); cfns = 'Annual';                end
if fns(2); cfns = [cfns ' & semi-annual']; end
if fns(3); cfns = [cfns ' & Time'];        end
if fns(4); cfns = [cfns ' & Time^2'];      end

if bogus>1 & (~fns(1) & ~fns(2))
   disp('LODRIVE Warning: Temporal bogusing disabled because not fitting')
   disp('temporal fns at this depth.')
   bogus = 1;
end   
   


% -------- Prepare the data ---------------------------------------

xd = xd(:);  yd = yd(:); dd = dd(:); tims = tims(:);

bdep = []; gdep = [];

% Prepare for bathymetry-influenced weighting system
% First, get any missing depths.

if TAR
   bdep = -get_bath(xd,yd);
   minbdep = 25;
   maxbdep = 2000;
   bdep = max([repmat(minbdep,size(bdep)) bdep]')';
   bdep = min([repmat(maxbdep,size(bdep)) bdep]')';
else
   bdep = zeros(size(xd));
end

lowrap = (any(xd<30) & any(xd>330));

if ~isempty(fns) && any(fns==1)
   doy = time2doy(tims);
   if ~(fns(3) || fns(4))
      tims = doy;
   end
else
   doy = zeros(size(xd));
   tims = doy;
end

% Use the next available name for output file
outputfl = 1;
if length(fnm)==1 && fnm==0
   if nargout==0
      % No output file nor arguments, which is dumb, so overrule and give
      % output file the default name.
      fnm = [];
   else
      outputfl = 0;
   end
end
if isempty(fnm)
   [tmp,fnm] = next_name(['lmap_' num2str(depm)]);
end

% Now get ocean depths at mapping grid points

if TAR
   gdep = -get_bath(xg,yg);
   if smth_dep
      % Smooth bathymetry (to smooth effects of TAR) using 9 point mean, triple
      % weight to centre point.
      [m,n] = size(xg);
      b = [gdep(:,1) gdep gdep(:,n)];
      b = [b(1,:); b; b(m,:)];
      y1=1:m; y2=2:m+1; y3=3:m+2;
      x1=1:n; x2=2:n+1; x3=3:n+2;
      gdep = (b(y1,x1)+b(y1,x2)+b(y1,x3)+b(y2,x1)+b(y2,x3)+b(y3,x1)+...
	      b(y3,x2)+b(y3,x3) + 3*b(y2,x2))/11;
   end
end

% Record configuration
config.depm = depm;
config.q1 = q1;
config.q2 = q2;
config.rmax = rmax;
if var_rmin
   config.rmin = [];
else
   config.rmin = rmin;
end
config.var_rmin = var_rmin;
if var_xw
   config.xw = [];
else
   config.xw = xw;
end
config.var_xw = var_xw;
config.bthy = TAR;
config.bogus = bogus;
config.fns = fns;
config.bar = bar;
if bar
   sbar = ' BAR system ';
else
   sbar = ' ';
end

config.m2coast = m2coast;
config.smw = smw;
config.smth_dep = smth_dep;
config.interpfill = interpfill;
config.lowrap = lowrap;
config.date = date;
config.constraint = si;
config.auto_ext_chnk = auto_ext;
config.noextrap = noextrap;
config.rQCn = rQCn;
config.wgtest = wgtest;
config.comment = cmt;
if ~isempty(damp)
   if ~isempty(damp{1}); config.damp_fn1 = damp{1}; end
   if ~isempty(damp{2}); config.damp_fn2 = damp{2}; end
   if ~isempty(damp{3}); config.damp_fn3 = damp{3}; end
   if ~isempty(damp{4}); config.damp_fn4 = damp{4}; end
end	 

if auto_ext  
   xm = []; ym = [];
   chnkstr = ' map_chunks auto-extend, '; 
else
   chnkstr = [' map_chunks [ym xm]=' num2str([ym xm])];       
end
if var_rmin
   rminstr = 'variable';
else
   rminstr = num2str(rmin);
end
if var_xw
   xwstr = ' xw=variable';
else
   xwstr = [' xweight=' num2str(xw)];
end
if TAR
   btstr = [', TAR ' num2str(TAR) ' '];
   if smth_dep
      btstr = [btstr '- pre-smoothed, '];
   end
else
   btstr = ', NO TAR ';
end
switch bogus
  case -4
    bgstr = ' spatial & temporal bogus to zero ';
  case -3
    bgstr = ' Bogus zero field strongly applied ';
  case -2
    bgstr = ' Bogus zero field applied ';
  case -1
    bgstr = ' Bogus zero field weakly applied ';
  case 0
    bgstr = ' Bogussing disabled ';       
  case 1
    bgstr = ' Spatial Bogus mean weight=exp(-(1.2*r/r(Q/10)).^3) ';
  case 2
    bgstr = ' Spatial & temporal bogussing ';
  case 3
    bgstr = ' Spatial & strong temporal bogussing ';
  case 4
    bgstr = ' Weak spatial & strong temporal bogussing ';
end       

igd = find(~isnan(dd));
dd = dd(igd);
xd = xd(igd);
yd = yd(igd);
tims = tims(igd);
doy = doy(igd);
if ~isempty(dw)
   dw = dw(igd);
end
if ~isempty(bdep)
   bdep = bdep(igd);
end

details=['Depm ' num2str(depm) ' ' num2str(length(igd),10) ' casts, '...
	 cfns ' q ' num2str([q1 q2]) ', rmin=' rminstr...
	 ' rmax=' num2str(rmax) ' varsort ' chnkstr ...
	 xwstr ', w=(1-r^3)^3' btstr bgstr sbar ', LOESS_MAP ' date];

if ~isempty(cmt)
   details = [details ' Comment: "' cmt '"'];
end


% Prepare to save outputs, and do mapping, and optionally fill gaps

an = []; sa = [];
if outputfl
   scmd=['save ' fnm ' config details q1 q2 rmax mn rq nq'];
   if fns(1);  scmd = [scmd ' an'];    end
   if fns(2);  scmd = [scmd ' sa'];    end
   if fns(3);  scmd = [scmd ' tco'];   end
   if fns(4);  scmd = [scmd ' t2co'];  end
   if smw;     scmd = [scmd ' swgt'];  end
end

[mn,rq,nq,an,sa,swgt,tco,t2co] = map_chunks(xd,yd,tims,dd,dw,bdep,...
					    xg,yg,gdep,config,xm,ym,rmin,xw); 


nodat = find(nq<=0);   
isdat = find(nq>0 & ~isnan(nq));
if interpfill & ~isempty(nodat)
   [mn,an,sa] = gapfiller(xg,yg,nodat,isdat,15,[],mn,an,sa);
   details = [details ' Interp ' num2str(length(nodat)) ' gaps'];
end
disp(details)

if outputfl
   if all(xg(1,:)==xg(end,:)) & all(yg(:,1)==yg(:,end)) 
      % If plaid, grid is completely specifed by 2 vectors
      Xg = xg(1,:);
      Yg = yg(:,1)';
      scmd = [scmd ' Xg Yg'];
   else
      scmd = [scmd ' xg yg'];
   end     
end


% -------------- Get stats fields, for QC screening or not -------------------
if scrn>0
   % In the following, we assume that nans are present in the same places in
   % the mean and the temporal fields. If not, need to go back to the slower
   % old method of individual calls to gapfiller for each field.

   jj = find(isnan(mn));
   ii = find(~isnan(mn));
   [mne,ane,sae] = gapfiller(xg,yg,jj,ii,[],[],mn,an,sa);
   res = dd - atdaypos(yd,xd,doy,xg,yg,mne);
   if ~fns(1) & ~fns(2)
      ii = find(~isnan(res));
      [sd,rmsmr,bcnt] = res_chunks(xd(ii),yd(ii),dd(ii),res(ii),xg,yg,[]);
      scmd = [scmd ' sd rmsmr'];
   else
      if ~fns(2)
	 res2 = dd - atdaypos(yd,xd,doy,xg,yg,mne,ane);
      else
	 res2 = dd - atdaypos(yd,xd,doy,xg,yg,mne,ane,sae);
      end

      ii = find(~isnan(res) & ~isnan(res2));
      [sd,rmsmr,bcnt,rmsr] = res_chunks(xd(ii),yd(ii),dd(ii),res(ii),xg,yg,res2(ii));
      res = res2;
      if fns(1) & any(abs(ane)>3.*sd)
	 disp(['NOTE: abs(an) exceeds 3xSD: ' num2str(sum(abs(ane)>3.*sd))]);
      end
      if fns(2) & any(abs(sae)>3.*sd)
	 disp(['NOTE: abs(sa) exceeds 3xSD: ' num2str(sum(abs(sae)>3.*sd))]);
      end      
      scmd = [scmd ' sd rmsmr rmsr'];
   end
   
   
   disp(['Data    : ' num2str([nanmin(dd) nanmax(dd) ...
		    nanmean(dd) nanstd(dd)])]);
   disp(['RMS_RES : ' num2str([nanmin(rmsr(:)) nanmax(rmsr(:)) ...
		    nanmean(rmsr(:)) nanstd(rmsr(:))])]);
   disp(['Minimum allowed rms_res: ' num2str(rmsmin) '. Reject > ' ...
	 num2str(scrthr) '*rms_res']);
   disp(['ABS     : ' num2str([nanmin(abs(res)) nanmax(abs(res)) ...
		    nanmean(abs(res)) nanstd(abs(res))])]);   
end


% ---------- Use stats fields to do QC screening -----------------

if scrn==1 && any(~isnan(rmsr(:)))
   jj = find(rmsr<rmsmin);
   if ~isempty(jj)
      rmsr(jj) = repmat(rmsmin,size(jj));
      disp(['Min rms_res asserted ' num2str(length(jj),6) ' times.']);
   end
   
   rmsr = spreadNcell(rmsr,2,0);
   
   rmsint = interp2(xg,yg,rmsr,xd,yd,'*linear');
   
   rejidx = find(abs(res./rmsint) > scrthr );

   rejpc = length(rejidx)*100/sum(~isnan(res));
   minrms = nanmin(rmsint(:));
   scrrep = ['Min rmsint= ' num2str(minrms) ',  ' ...
	     num2str(rejpc) ' % > ' num2str(scrthr) ' rms(residuals)'];
   disp(scrrep);   

   scmd = [scmd ' rejidx'];
end

if outputfl
   eval(scmd);
end
   
if ~isempty(ppname)
   eval(['!rm ' ppname '*.mat']);
end

return

%-------------------------------------------------------------------------

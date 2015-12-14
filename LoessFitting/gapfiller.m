% GAPFILLER   Interpolate to fill gaps in matrix.
%
% Works for complex as well as real values.  
% Note that this is much quicker than fillgap because works on
% small chunks. It will do poorly (especially at small bsz, such as the
% default 5) at filling large gaps but is actually much better (more local)
% for infilling onto a higher resolution grid.
%
% nodat and isdat need not comprise a complete set - there may be a 3rd class 
% of points (don't want to fill them, but no good data there).
%
% bsz:  bin size in number of grid intervals (also determines min number of points) 
%       Bin expands until have (bsz^2)/3 good points
%
% minr: Supply if only want to fill where have at least one datum within a
%       min radius (in degrees)
%
% USAGE: [fil1,fil2,...] = gapfiller(xx,yy,nodat,isdat,bsz,minr,dat1,dat2,...);

function varargout = gapfiller(xx,yy,nodat,isdat,bsz,minr,varargin)

% MODS: 23/10/07  minr didn't work before - fixed

if nargin<6
   minr = [];
end
if nargin<5 | isempty(bsz)
   bsz = 5;
end
mind = floor((bsz^2)/3);

nv = nargin-6;
if nargout~=nv
   error('Need same number of input and output data variables')
end

dud = [];
vars = 1:nv;
for ii = vars
   varargout{ii} = varargin{ii};
   if isempty(varargin{ii}); dud = [dud ii]; end
end
if ~isempty(dud)
   vars(dud) = [];
end

if isempty(nodat) | length(isdat)<mind
   return
end

[ny,nx] = size(xx);
[cx,cy] = meshgrid(1:bsz:nx,1:bsz:ny);
mdim = max([nx ny]);

idx = zeros(ny,nx);
idx(isdat) = isdat;
idx(nodat) = -nodat;

if ~isempty(minr)
   % Pre-square minr because we test against x^2+y^2, not sqrt(x^2+y^2)
   minr = minr^2;
end

for kk = 1:prod(size(cx))
   gd = [];
   ext = 0;
   while length(gd) < mind  &  ext < mdim
      ext = ext+bsz;
      xch = max([1 cx(kk)-ext]):min([cx(kk)+(bsz-1)+ext nx]); 
      ych = max([1 cy(kk)-ext]):min([cy(kk)+(bsz-1)+ext ny]); 

      cidx = idx(ych,xch);
      gd = cidx(find(cidx>0));
   end
   
   xin = cx(kk):min([cx(kk)+(bsz-1) nx]); 
   yin = cy(kk):min([cy(kk)+(bsz-1) ny]); 

   cidx = idx(yin,xin);
   ii = -cidx(find(cidx<0));

   for jj = ii(:)'
      x = xx(gd)-xx(jj);
      y = yy(gd)-yy(jj);

      % 4th power - slower but gives more weight to local data 
      %   r = x.^2.*x.^2 + y.^2.*y.^2;     % x^4+y^4, but vastly faster as coded

      r = x.^2 + y.^2;
      ok = (isempty(minr) || min(r)<minr);
      r = r.^2;
      rsum = sum(1./r);
      if ok
	 for ll = vars
	    varargout{ll}(jj) = sum(varargin{ll}(gd)./r)./rsum;
	 end
      end
   end
end

%for jj = nodat(:)'
%   x = abs(xx(gd)-xx(jj));
%   y = abs(yy(gd)-yy(jj));
%
%   % 4th power - slower but gives more weight to local data 
%   r = 1./(x.^2.*x.^2 + y.^2.*y.^2);     % x^4+y^4, but vastly faster as coded
%
%   % 3rd power - reverts to global mean more quickly as move into data void.
%   %  r = 1./(x.*x.^2 + y.*y.^2);     % x^3+y^3, vastly faster as coded
%
%   fil1(jj) = sum(dat1(gd).*(r))./sum(r);
%   if nargout>1
%      fil2(jj) = sum(dat2(gd).*(r))./sum(r);
%   end
%   if nargout>2
%      fil3(jj) = sum(dat3(gd).*(r))./sum(r);
%   end
%end

%===========================================================================

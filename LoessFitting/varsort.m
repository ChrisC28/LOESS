% VARSORT  -  MEX function to "sort" for loess mapping. Used where 
% require the q smallest values, but they don't have to be in order. 
% Varsort estimates the value of the q-th value, without sorting, so that
% the matlab routine can then use:   Q = find(r<rq)
% Varsort differs from earlier "binsort" in that the number required 
% decreases as value (radius) increases. It is also about 40% slower.
%   
%   Call:
%           [rq,nq] = varsort(r,q1,q2,maxr,minr)
%  where:
%      r  - vector of values (radii) to "sort"
%      q1 - required number of values to select at minimum radius
%      q2 - required number of values to select at maximum radius
%      maxr - maximum radius
%      minr - minimum radius
%
%   This m-file for VARSORT is also available (it is fashioned with minimal
%   change from the C code, to enable careful testing of the C code.)

function [rad,nsel] = varsort(orig,q1,q2,maxr,minr)

nval = length(orig);
NCAT = 100;
SLACK = 3;
   
% Find max value - we actually do this much faster than Matlab! 
% Return immediately if have degenerate situation of all points being 
% inside min r.
% Increase rmax by a fraction so furthest point lies within the last cell

rmax = max(orig);

if rmax<minr
   nsel = nval;
   rad = rmax*1.0001;
   return;
elseif rmax>maxr
   rmax = maxr;
end
   
rmax = rmax*1.001;
rrng = rmax-minr;

% Zero bin count and construct number_required_vs_distance vector

rtmp = sqrt(q1-q2)/(NCAT-1);
for kk=1:NCAT 
   cnt(kk) = 0;   
   nreq(kk) = q1 - (rtmp*(kk-1)).^2;
end

% Count number in each bin

cnt0 = 0;
for ii=1:nval
   if orig(ii)<minr
      cnt0=cnt0+1;
   elseif orig(ii)<rmax
      jj = 1+floor(NCAT*(orig(ii)-minr)/rrng);
      cnt(jj) = cnt(jj)+1;
   end
end

% If already have required number in min radius, return with that radius

if cnt0>=q1
   nsel = cnt0;
   rad = minr;
   return;
end
   
% Step through until aggregate counts > number required

nsel = cnt0+cnt(1);
i1 = 1;   
while nsel<(nreq(i1)-SLACK) & i1<NCAT
   i1=i1+1;
   nsel = nsel+cnt(i1);
end
   
% If necessary (no. selected not close enough to no. required), do a second
% pass, subdividing the bin containing the nreq-th element

if nsel<=(nreq(i1)+SLACK)
   rad = minr+((i1+1)*(rrng/NCAT));
else
   nsel = nsel-cnt(i1);
   subsiz = rrng/NCAT;
   cmin = minr + ((i1-1)*subsiz);
   cmax = cmin+subsiz;
      
   for kk=1:NCAT; cnt(kk) = 0; end

   for ii=1:nval
      if orig(ii)>cmin & orig(ii)<cmax
	 jj = 1+floor(NCAT*(orig(ii)-cmin)/subsiz); 
	 cnt(jj) = cnt(jj)+1;
      end
   end

   i2 = 0; 
   while nsel<(nreq(i1)-SLACK) & i2<NCAT
      i2=i2+1;
      nsel = nsel+cnt(i2);
   end

   rtmp = cmin + i2*(subsiz/NCAT);
   
   % Interp removed because actually bad where have data clumps!
   %% Perhaps an unnecessary correction, but lin interp a value within the
   %% last category 
   %
   %if nsel>(nreq(i1)+SLACK) & cnt(i2)>1
   %   rtmp = rtmp - (subsiz/NCAT)*((nsel-nreq(i1))/cnt(i2));
   %   nsel = nreq(i1);
   %end
      
   rad = rtmp;
end

return;


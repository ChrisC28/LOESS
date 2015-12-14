function [siz]=binsort(orig,nreq)

% BINSORT Find smallest nreq values in a vector. NOT USED now because
%          binsort.c is more efficient.
% INPUT: orig - vector of values (ie radii in km)
%        nreq - number of values required
%
% OLD Call sequence:
%    if first==1
%       first=0;
%       rqj = max(r)/2;
%    end
%    [Q,rqj] = binsort(r,q,rqj);
%
% USAGE:   rqj = binsort(r,q);

% Jeff Dunn  11/4/96
% May be sufficient, or in some cases best, just use pure binary-sort.
% For now, refine by using:
% a)  size limit determined from max length by a ratio of num_req/tot_num.
% b)  input last radius as first guess

slack = 2;

% NOTE:   Premature RETURN for degenerate case
if length(orig) < nreq
 siz = max(orig);
 return
end

piece = orig;
ntot = 0;  mn = 0;
mx = max(piece);
inc = mx;
mininc = mx/10000;     % A threshold to stop search if inc gets too small
nmore = nreq;
siz = mean(piece);

while nmore > slack
 sel = find(piece < siz);
 nsel = length(sel);
 if inc<mininc
   if nsel<nmore
     siz = siz+inc/2;
   else
     siz = siz-inc/2;
   end
   nmore = 0;
 elseif abs(nsel-nmore)<=slack
   nmore = 0;
 else
   if nsel<nmore
     ntot = ntot+nsel;
     nmore = nreq-ntot;
     piece(sel) = [];
     mn = min(piece);
   else
     piece = piece(sel);
     mx = siz;
   end
   inc = (mx-mn)*nmore/length(piece);
   siz = mn + inc;
 end

end

% indx = find(orig < siz);

%if abs(length(indx)-nreq)>6
%  disp(['### ' num2str(length(indx)) ' ' num2str(inc,7) ' ' num2str(nsel)]);
%end

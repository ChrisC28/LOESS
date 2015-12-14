% Find the last filename in the series
function [flast,fnext] = next_name(fbase,fpath)

if nargin<2 | isempty(fpath)
   fpath = '';
end

flast = [];
ich = 64;
fnext = fbase;
while exist([fpath fnext '.mat'],'file')
   flast = fnext;
   ich = ich+1;
   fnext = [fbase char(ich)];
end

return

% MAXJD: On machines where the mex maxjd is invoked, it is much more efficient
%        than the Matlab 5 builtin function. On other machines, maxjd is just
%        a wrapper for the Matlab "max".
%
% NOTE:  Returns a NaN for any column containing a NaN.
%
% USAGE: maxA = maxjd(A)     OR    [maxA,indx] = maxjd(A)

% This code is only as a fallback on machines on which I cannot compile maxjd.c

function [ma,indx] = maxjd(a)

if nargout==1
  ma = max(a);
else
  [ma,indx] = max(a);
end

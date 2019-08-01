function [empty,a] = emptyintersect(a,b)
%EMPTYINTERSECT    Compute intersection and check for empty components
%
%This is basically a dummy routine for floating-point input.
%
%   [empty,c] = emptyintersect(a,b)
%
%Result c is that of intersect(a,b), and 
%  empty(i) = 1     a(i)~=b(i)
%             0     a(i)==b(i)
%             NaN   at least one of a(i) and b(i) is NaN
%

% written  04/23/12     S.M. Rump  (Thanks to Kolumb�n S�ndor for pointing
%                                     to that missing function)
% modified 07/16/15     S.M. Rump  comment
%

  empty = ( a~=b );         % Note x~=NaN is always true
  a(find(empty)) = NaN;
  
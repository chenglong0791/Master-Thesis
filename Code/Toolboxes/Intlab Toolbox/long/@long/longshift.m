function A = longshift(A,r)
%LONGSHIFT    Shifts long number A by r bits
%
%   A = longshift(A,r)
%
%The same as A * long(2)^r
%

% written  01/10/04     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
% modified 04/23/14     S.M. Rump  set/getappdata replaced by global
% modified 07/30/16     S.M. Rump  rounding check by getround for Matlab 2016b 
%

  global INTLAB_CONST

  rndold = getround;
  if rndold
    setround(0)
  end

  INTLAB_LONG_LOGBETA = INTLAB_CONST.LONG_LOGBETA;

  q = floor(r/INTLAB_LONG_LOGBETA);
  A = struct(A);
  A.exponent = A.exponent + q;
  rem = r - q*INTLAB_LONG_LOGBETA;
  A.mantissa = pow2(A.mantissa,rem);
  A = normalize(A);
  A = normalizefirst(A);
  A = class(A,'long');
    
  if rndold
    setround(rndold)
  end

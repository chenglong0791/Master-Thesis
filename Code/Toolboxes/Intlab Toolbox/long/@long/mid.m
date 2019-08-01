function C = mid(A)
%MID          Midpoint of long (only with error term)
%
%  C = mid(A)
%

% written  11/06/99     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
% modified 04/23/14     S.M. Rump  set/getappdata replaced by global
%

  global INTLAB_CONST

  INTLAB_LONG_ERROR = INTLAB_CONST.LONG_ERROR;
  
  if ~INTLAB_LONG_ERROR
    error('mid called for long without error term')
  end

  C = A;
  C.error.mant = 0;
  C.error.exp = 0;

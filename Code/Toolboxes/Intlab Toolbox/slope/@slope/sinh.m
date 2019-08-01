function u = sinh(a)
%SINH         Slope hyperbolic sine sinh(a)
%

% written  12/06/98     S.M. Rump
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

  INTLAB_SLOPE = INTLAB_CONST.SLOPE;

  u = a;

  u.r = sinh(a.r);
  indexc = 1:INTLAB_SLOPE.NUMVAR;
  indexr = 2:INTLAB_SLOPE.NUMVAR+1;
  Xxs = hull(a.r(:,indexc),a.r(:,indexr));
  Index = 1:size(a.r.inf,1);

  index = all( a.r.sup<=0 , 2);
  if any(index)
    aindex.r = a.r(index,:);
    aindex.s = a.s(index,:);
    u.s(index,:) = slopeconvexconcave('sinh','cosh(%)',aindex,0);
    Index(index) = 0;
  end

  index = all( a.r.inf>=0 , 2);
  if any(index)
    aindex.r = a.r(index,:);
    aindex.s = a.s(index,:);
    u.s(index,:) = slopeconvexconcave('sinh','cosh(%)',aindex,1);
    Index(index) = 0;
  end

  if any(Index)
    Index( Index==0 ) = [];
    u.s(Index,:) = a.s(Index,:) .* cosh(Xxs(Index));
  end
  
  if rndold
    setround(rndold)
  end

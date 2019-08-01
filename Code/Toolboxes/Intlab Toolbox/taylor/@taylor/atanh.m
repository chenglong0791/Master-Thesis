function r = atanh(a)
%ATANH        Taylor inverse hyperbolic tangent  atanh(a)
%
%Thanks to George Corliss for providing the Taylor expansion
%

% written  06/03/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
% modified 04/23/14     S.M. Rump  set/getappdata replaced by global
% modified 07/30/16     S.M. Rump  rounding check by getround for Matlab 2016b 
%

  global INTLAB_CONST
  
  rndold = getround;
  if rndold
    setround(0)
  end

  K = INTLAB_CONST.TAYLOR_ORDER;
  
  r = a;
  ct = a.t;
  N = size(a.t,2);
  r.t(1,:) = atanh(a.t(1,:));
  ct1 = 1-a.t(1,:).^2;     % 1-a^2
  r.t(2,:) = a.t(2,:) ./ ct1 ;
  for j=2:K
    ct(j,:) = sum( a.t(1:j,:).*a.t(j:-1:1,:) , 1 );
    r.t(j+1,:) = ( j*a.t(j+1,:) + sum( repmat((1:j-1)',1,N).*r.t(2:j,:).*ct(j:-1:2,:) , 1 ) ) ./ ( j*ct1 );
  end

  if rndold
    setround(rndold)
  end

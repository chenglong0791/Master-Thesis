function c = midrad(a,r)
%MIDRAD       Initialization of interval by midpoint and radius
%  computed such that <a,r> is enclosed in interval  c
%
%   c = midrad(a,r)
%

% written  10/16/98     S.M. Rump
% modified 06/24/99     S.M. Rump  check sparsity, multi-dimensional arrays
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 11/16/02     S.M. Rump  sparse radius
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  faster check for rounding to nearest
% modified 09/10/07     S.M. Rump  performance
% modified 10/14/08     S.M. Rump  huge midpoint array and 0 radius
% modified 10/14/08     S.M. Rump  transform [inf,inf] into [realmax,inf] or <0,inf>
% modified 02/28/12     S.M. Rump  check input is numeric (thanks to Joshua Ross, Adelaide)
% modified 07/11/12     S.M. Rump  typo
% modified 05/15/14     S.M. Rump  code optimization
% modified 12/09/15     S.M. Rump  prod(size) to numel
% modified 07/30/16     S.M. Rump  rounding check by getround for Matlab 2016b 
%

  if ~isnumeric(a) || ~isnumeric(r)
    error('both bounds for infsup must be numeric')
  end

  index = any( r<0 );                   % take care of huge matrices
  if ( ~isreal(r) ) || any(index(:))
    error('invalid radius in call of midrad')
  end

  if ( numel(r)~=1) && ( issparse(a)~=issparse(r) )
    error('midpoint and radius must both be or none of them be sparse')
  end

  if isreal(a)                     % real interval (inf/sup representation)
    rndold = getround;
    if issparse(a) && isequal(r,0)   % already proper intervals
      r = sparse(0);
    end
    setround(-1)
    cinf = a - r;
    setround(1)
    csup = a + r; 
    setround(rndold)                % set rounding to previous value
    c = intval(cinf,csup,'infsup');
  else                             % complex interval (mid/rad representation)
    c = intval(a,r,'midrad');
  end

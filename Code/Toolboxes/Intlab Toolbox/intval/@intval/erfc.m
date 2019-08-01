function y = erfc(x)
%ERFC         Implements  erfc(x)  for real intervals
%
%   y = erfc(x)
%
%interval standard function implementation
%

% written  05/30/13     S.M. Rump
% modified 07/30/16     S.M. Rump  rounding check by getround for Matlab 2016b 
%

  if x.complex
    error('Error function only for real arguments')
  end
  
  if issparse(x)
    index = ( x==0 );
    if any(index(:))                    % treat zero indices
      y = intval(ones(size(x)));
      index = ~index;
      %VVVV  y(index) = erfc(full(x(index)));
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,erfc(full(subsref(x,s))));
      %AAAA  Matlab bug fix
    else
      y = erfc(full(x));
    end
    return
  end

  rndold = getround;
  if rndold
    setround(0)
  end
  
  y = x;
  
  if all( x.inf(:)==x.sup(:) );         % thin input
    
    [y.inf,y.sup] = erfc_rnd(x.inf(:),[]);
    
  else                                  % thick input
    
    y.inf = erfc_rnd(x.sup(:),-1);
    y.sup = erfc_rnd(x.inf(:),1);
    
  end
  
  y.inf = reshape(y.inf,size(x.inf));
  y.sup = reshape(y.sup,size(x.sup));
  
  setround(rndold)  

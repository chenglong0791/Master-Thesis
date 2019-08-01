function setting = intlabsetting(setting,see)
%SETTING      Current setting of INTLAB control variables
%
%For display of current setting, call
%
%   intlabsetting(see)
%
%Parameter see is optional; default is 1 for displaying results. For storing
%and restoring current setting use
%
%   currentsetting = intlabsetting  [or  currentsetting = intlabsetting(see)]
%   intlabsetting(currentsetting)   [or  intlabsetting(currentsetting,see)  ]
%
%For restoring the default as specified in startintlab use
%
%   intlabsetting('default')
%

% written  04/04/04     S.M. Rump
% modified 02/12/06     S.M. Rump  new settings added
% modified 09/10/07     S.M. Rump  approximate std fcts removed
% modified 10/23/07     S.M. Rump  new setting corrected
% modified 11/24/07     S.M. Rump  isstr replaced by ischar
% modified 10/18/08     S.M. Rump  out-of-range flag
% modified 09/20/10     S.M. Rump  global variables declaration first
% modified 08/24/12     S.M. Rump  STAGE superflous due to redesign of verifylss
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/03/12     S.M. Rump  SparseArrayDeriv for gradients, hessian and slope removed
% modified 04/23/14     S.M. Rump  set/getappdata replaced by global
% modified 04/30/14     S.M. Rump  fl- and affari package
% modified 05/15/14     S.M. Rump  code optimization
%

  global INTLAB_CONST
  
  INTLAB_SETTING = INTLAB_CONST.SETTING;
  % flag: ignore input out of range; for internal use only
  INTLAB_CONST.RealStdFctsExcptnIgnore = 0; % default setting: inactive

  if ( nargin>=1 ) && ischar(setting)    % new setting
    if isequal(setting,'default')        % restore default setting      
      setting = INTLAB_SETTING;
    end
    if nargin==1
      see = 1;
    end
    for i=1:size(setting,1)
      j = strfind(setting(i,:),')');
      eval([ setting(i,1:j-1) ',0' setting(i,j:end) ';' ]);
    end
    if see
      intlabsetting
    end
    return
  end
  
  if nargin==0
    see = 1;            % display and possibly output of current setting
  else
    see = setting;
  end
  
  digits = displaywidth; setting = ['displaywidth(' int2str(digits) ');'];
  
  s = intvalinit('Display'); setting = char( setting , ['intvalinit(''' s ''');'] );
  s = intvalinit('RealStdFctsExcptn'); setting = char( setting , ['intvalinit(''' s ''');'] );
  s = intvalinit('Residual'); setting = char( setting , ['intvalinit(''' s ''');'] );
  s = intvalinit('IVmult'); setting = char( setting , ['intvalinit(''' s ''');'] );
  s = intvalinit('RealComplexAssign'); setting = char( setting , ['intvalinit(''' s ''');'] );
  s = intvalinit('ComplexInfSupAssign'); setting = char( setting , ['intvalinit(''' s ''');'] );
  
  s = affariinit('Display'); setting = char( setting , ['affariinit(''' s ''');'] );
  s = affariinit('RoundingErrors'); setting = char( setting , ['affariinit(''' s ''');'] );
  s = affariinit('ApproxMode'); setting = char( setting , ['affariinit(''' s ''');'] );

  s = flinit('Display'); setting = char( setting , ['flinit(''' s ''');'] );
  s = flinit('AccumulateMode'); setting = char( setting , ['flinit(''' s ''');'] );
    
  s = polynominit('DisplayUPoly'); setting = char( setting , ['polynominit(''' s ''');'] );
  s = polynominit('EvaluateUPoly'); setting = char( setting , ['polynominit(''' s ''');'] );
  s = polynominit('EvaluateMPoly'); setting = char( setting , ['polynominit(''' s ''');'] );
  s = polynominit('AccessVariable'); setting = char( setting , ['polynominit(''' s ''');'] );
  
  s = longinit('ErrorTerm'); setting = char( setting , ['longinit(''' s ''');'] );
  
  if see
    for i=1:size(setting,1)
      eval(setting(i,:))
    end
  end
  
  if nargout==0
    clear setting
  end
  
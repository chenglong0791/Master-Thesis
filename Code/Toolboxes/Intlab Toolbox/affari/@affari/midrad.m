function midrad(c)
%MIDRAD       Display affaris in midrad notation
%
%   midrad(c)
%

% written  04/04/14     S.M. Rump
% modified 04/23/14     S.M. Rump  set/getappdata replaced by global
% modified 05/17/14     S.M. Rump  code optimization
%

  global INTLAB_CONST

  name = inputname(1);
  if isempty(name)                    % happens for midrad(affari(random))
    name = 'ans';
  end

  INTLAB_INTVAL_DISPLAY = INTLAB_CONST.INTVAL_DISPLAY;
  INTLAB_CONST.INTVAL_DISPLAY = 'DisplayMidrad';
  display(c,name)
  INTLAB_CONST.INTVAL_DISPLAY = INTLAB_INTVAL_DISPLAY;
      
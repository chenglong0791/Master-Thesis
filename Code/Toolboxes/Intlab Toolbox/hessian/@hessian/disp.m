function disp(x)
%DISP         Display function for pop-up windows in debugger
%

% written  04/26/13     S.M. Rump
% modified 10/17/16     S.M. Rump  restricted output for popup windows
%

  global INTLAB_CONST
  
  if numels(x)<INTLAB_CONST.DISPLAY_POPUP
    display(x);
  else
    disp('Display suppressed.')
  end
  
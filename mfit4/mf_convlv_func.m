function [y, name, pnames, pin]=mf_convlv_func(x,p, flag)
% function to convolute an other function with some data stored in mf_convlv window

fitfun=get(findobj('tag','mf_FitFuncFile'),'string');

if nargin==2
  y = feval(fitfun, x, p);
  % check if 1D convolution must be carried out
  hmf_convlv = findobj('Tag', 'mf_convlv');
  if ~isempty(hmf_convlv)
    y = mf_convlv('convlv', y);
  end
else
  [y, name, pnames, pin]=feval(fitfun, x, p, flag);
  name = [ name ' (convoluted)' ];
end
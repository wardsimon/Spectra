function mf_newfn(cmd)
%
% MFIT  function mf_newfn
%     Load new fitting function
%     MZ 29.11.94
%
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

if nargin == 0, cmd=''; end
if isempty(cmd) | strcmp(cmd,'noprompt')

%---------- Get fit function details ---------------
fundir=get(findobj('tag','mf_FitFuncDir'),'string');
funfile=get(findobj('tag','mf_FitFuncFile'),'string');
funname=get(findobj('tag','mf_FitFuncName'),'string');

else

funfile = cmd;
funname = cmd;
fundir = '';
set(findobj('tag','mf_FitFuncFile'),'string',cmd);
set(findobj('tag','mf_FitFuncName'),'string',cmd);
set(findobj('tag','mf_FitFuncDir'),'string','');

end

[a,b,pfix] = mf_rpars;
if ~isempty(hmf_pars) & hmf_pars
  set(findobj('tag','ChiSq'),'Userdata', [ a(:) b(:) pfix(:) ]);
end

%----------- Call fit func to initialize ------------------
% Call fit function with flag=1
pin = zeros(50,1);
[y, name, pnames, pin]=feval(funfile,0,pin,1);
npars=length(pin);

%----------- Create parameters window --------------------
prev_pars =hmf_pars;
hmf_pars=mf_pwin(funfile, pnames, pin, 1);
if (length(pin) == length(pfix)) & ~isempty(prev_pars) & prev_pars
  mf_upars([],[],pfix);
end

disp(['*Loaded fitting function: ' funname ' ('  funfile ').' ] )
mf_uplot('fit');

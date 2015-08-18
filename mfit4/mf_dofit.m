function mf_dofit
%
% MFIT function  mf_dofit
%   Do fit
%     MZ 29.11.94 EF 04.97
%
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

[p, dp, fixed]=mf_rpars;
if isempty(p)
  return;
elseif ~isempty(hmf_pars) & hmf_pars
  set(findobj('tag','ChiSq'),'Userdata', [ p(:) dp(:) fixed(:) ]);
end
npars=length(p);


%-------------- Extract fit function name and dir ------------------
fitfun=get(findobj('tag','mf_FitFuncFile'),'string');
fitrt =get(findobj('tag','mf_FitRoutineFile'),'string');
if isempty(fitrt)
  g=get(findobj('tag','mf_FitRoutineMenu'),'children');
  eval(get(g(end-1),'callback'))
  fitrt =get(findobj('tag','mf_FitRoutineFile'),'string');
end


%-------------- Extract data from figure -----------------------
data=get(hmf_data,'userdata');

%----------- Create vectors of selected data points -----------------
index=data(:,4);
xs=data(:,1);
ys=data(:,2);
errs=data(:,3);
if sum(index)<=1
  disp('Can''t fit : no points !!');
  return;
end
if (sum(index)-sum(dp))<0
  mf_msg('Too few points selected - unconstrained - interpoling data');
  %return;
end
idx = find(isnan(ys) | isinf(ys) | (index == 0));
errs(idx) = 0;

% --------- Do the fit ------------------------

[y, name, pnames, pin]=feval(fitfun,[],p,1);

disp(sprintf('* Fit using function: %s (%s)\n', name, fitfun));

% handling convolution case
hmf_convlv = findobj('Tag', 'mf_convlv');
if ~isempty(hmf_convlv)
  UserData = get(hmf_convlv, 'UserData');
  if UserData.output.convlv & get(UserData.handles.checkbox, 'Value')
    fitfun = 'mf_convlv_func';
  end
end

mf_msg(sprintf('Now fitting... %i points, %i used',length(xs),length(xs)-length(idx)));
% check if 1D convolution must be carried out
hmf_convlv = findobj('Tag', 'mf_convlv');
if ~isempty(hmf_convlv)
  UserData = get(hmf_convlv, 'UserData');
  mf_msg(['Convoluting fit function with ' UserData.input.name ' (' UserData.input.type ')']);
  mf_convlv('check', xs);
end
[p,dp]=feval(fitrt,xs,ys,errs,p,~fixed,fitfun);

datadir = get(findobj('Tag','mf_DataDir'),'String');          % data directory
datafile = get(findobj('Tag','mf_DataFile'),'String');
time = clock;
disp(['* Data : ' datadir datafile ]);
fprintf(1,'* Fit results: %s %i:%i\n', date, time(4), time(5))
disp('---------------------------------------------')
disp('  Parameter       Value      Uncertainty')

for i=1:length(p)
   disp(sprintf('%12s  %12.4e  %12.4e', pnames(i,:), p(i), dp(i)))
end
disp(sprintf('---------------------------------------------\n'))

if strcmp(fitfun,'multifunc')
  modifp = [ get(findobj('Tag','mf_mfconst'),'String') ';' ];
  if length(modifp) > 1
    disp('Evaluating constrain field in MultiFunction')
    fprintf(1,'for parameter results : %s\n',modifp);
    eval(modifp);
  end
  set(findobj('Tag','mf_mfn'),'Userdata', [ p(:) dp(:) fixed(:) ]);
    % store associated starting parameters
end

mf_upars(p, dp);
mf_uplot('fit');
mf_stats;

h = findobj('tag','mf_ExecAfterFit');
if ~isempty(h)
  todo = get(h,'String');
  if ~isempty(todo)
    fprintf(1,'Executing : %s\n',todo);
    eval(todo);
  end
end

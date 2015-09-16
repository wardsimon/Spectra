function [x,y,err,selected,fit,p,dp,fixed]=fromfit(mode)
% MFIT data collect
% function [x,y,err,selected,fit,p,dp,fixed]=fromfit;
%     Get data from MFIT
% The other 'direct access' function is 'tomfit'
%     EF 03.07.97
%

curfig = get(0,'CurrentFigure');
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

if nargin == 0, mode=[]; end

if isempty(mode)
  mode = 'normal';
end

%---------- Attach data to userdata ------------------
x=[];y=[];err=[];p=[];dp=[];fit=[];selected=[];fixed=[];

if (isempty(hmf_ctrl) | ~hmf_ctrl) & strcmp(mode,'normal')
  disp('MFit is not running.')
  return
end

if strcmp(mode,'normal')
  disp('Collect data from MFIT windows');
end
if ~isempty(hmf_data) & hmf_data
  userdata=get(hmf_data,'userdata');
  if ~isempty(userdata)
    x=userdata(:,1);
    y=userdata(:,2);
    err=userdata(:,3);
    selected=userdata(:,4);
  end
end
if isempty(x) & strcmp(mode,'normal')
  disp('No user data');
end

hfit=findobj('Tag','mf_fitline');
if ~isempty(hfit)
  mf_msg('Getting fit results');
  yfit= get(hfit,'Ydata');
  xfit= get(hfit,'Xdata');
  if iscell(yfit)
    yfit = yfit{1};
    xfit = xfit{1};
  end
  fit=reshape(mf_interp(xfit,yfit,x),size(x));
elseif strcmp(mode,'normal')
  disp('No fit data');
end

if ~isempty(hmf_pars) & hmf_pars
  hm=get(hmf_pars,'userdata');
  p = [];
  dp = [];
  fixed = [];
  [n,c]=size(hm);
  for i=1:n
  h=hm(i,1);
    p = [ p str2num(get(h,'String')) ];
    h=hm(i,2);
    dp = [ dp str2num(get(h,'String')) ];
    h=hm(i,3);
    fixed = [ fixed get(h,'Value') ];
  end
elseif strcmp(mode,'normal')
  disp('No pars data');
end

figure(curfig);

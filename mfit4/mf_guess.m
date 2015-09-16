function mf_guess(ags)
%
% MFIT function mf_guess(ags)
% ags = autoguess : 1, user guess : 0
% 		Make guess for parameters
% 		MZ 29.11.94
% 
if nargin < 1
	ags = [];
end
hmf_data=findobj('tag','mf_DataWindow');
hmf_par=findobj('tag','mf_ParWindow');
if isempty(hmf_data)
   errordlg('No data window open','MFIT error;');
   return
else
	data=get(hmf_data,'userdata');
	xs=data(:,1);
	ys=data(:,2);
	index=data(:,4);
	if sum(index) <= 1
		disp('No data selected');
		return;
	end
	idx = find(~isnan(ys) & ~isinf(ys) & (index ~= 0));
	xs = xs(idx); ys = ys(idx);
end

fitfun=get(findobj('tag','mf_FitFuncFile'),'string');
if isempty(hmf_par)
	pin = zeros(50,1);
	[y, name, pnames, pin]=feval(fitfun,xs,pin,1);
	mf_pwin(name, pnames, pin);
end

%-------------- Extract fit function name and dir ---------------
p = mf_rpars;

%---------- Change to function dir, and call with flag=2 or autoguess --------

set(0,'CurrentFigure',hmf_data);
figure(hmf_data);
drawnow

[oldp, dp, fixed] = mf_rpars;
ifixed = find(fixed);

if (isempty(ags))
	ags = str2num(get(findobj('tag','mf_AutoGuess'),'string'));
end
if (ags)
	[pin,pmat] = mf_autogs; % keep fixed params
	set(findobj('Tag','mf_autoguess'),'userdata',pmat);
	pin = pin(:,1);
	[y, name, pnames, pin]=feval(fitfun,xs,pin,ys);
else
	refresh(hmf_data);
	figure(hmf_data);
	[y, name, pnames, pin]=feval(fitfun,xs,p,2);
end

if ~isempty(ifixed)
	pin(ifixed) = oldp(ifixed);	
end

%---------- Call pwin in case number of pars has changed --------
mf_pwin(name, pnames, pin);
sig=[];
mf_upars(pin, sig, fixed);
mf_uplot('fit');
